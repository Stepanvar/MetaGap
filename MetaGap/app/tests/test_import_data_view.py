"""Integration tests covering the VCF import view."""

from __future__ import annotations

from unittest import mock

from django.contrib.auth import get_user_model
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import TestCase
from django.urls import reverse

from ..models import AlleleFrequency, SampleGroup


class ImportDataViewIntegrationTests(TestCase):
    """Ensure VCF uploads produce related allele frequency records."""

    VCF_CONTENT = """##fileformat=VCFv4.2
##contig=<ID=1>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=CLNSIG,Number=1,Type=String,Description="Clinical significance">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##SAMPLE=<ID=GroupA,Description=Imported group,CustomKey=Value>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample001
1\t1234\trsTest\tA\tT\t99\tPASS\tAF=0.5;CLNSIG=Pathogenic\tGT:GQ\t0/1:99
"""

    VCF_WITH_METADATA = """##fileformat=VCFv4.2
##contig=<ID=1>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##SAMPLE=<ID=GroupB,Description=Detailed group>
##SAMPLE_GROUP=<Source_Lab=MetaLab,Contact_Email=lab@example.com,Contact_Phone=+123456789,Total_Samples=5,Inclusion=Adults,Exclusion=Under18>
##REFERENCE_GENOME_BUILD=<Build_Name=GRCh38,Build_Version=v1,Additional_Info={\"notes\":\"GRCh\"}>
##GENOME_COMPLEXITY=<Size=3.2Gb,Ploidy=Diploid,GC_Content=41%>
##SAMPLE_ORIGIN=<Tissue=Blood,Collection_Method=Venipuncture,Storage_Conditions=-80C,Time_Stored=6m>
##MATERIAL_TYPE=<Material_Type=DNA,Integrity_Number=9.5>
##LIBRARY_CONSTRUCTION=<Kit=HyperPrep,Fragmentation=Enzymatic,Adapter_Ligation_Efficiency=0.9,PCR_Cycles=8>
##SEQUENCING_PLATFORM=<Platform=Illumina NovaSeq,Instrument=NovaSeq 6000,Flow_Cell=FC123,Channel_Method=Two-Channel,Cluster_Density=1200k,QC_Software=Control,ReadLength=150>
##ONT_SEQ=<Instrument=PromethION,Flow_Cell=FLO-MIN106,Flow_Cell_Version=v1,Pore_Type=R9.4,Bias_Voltage=180mV>
##PACBIO_SEQ=<Instrument=Sequel II,Flow_Cell=SMRT Cell 8M,SMRT_Cell_Type=1M,ZMW_Density=High>
##IONTORRENT_SEQ=<Instrument=Ion S5,Flow_Cell=IonChip,Chip_Type=540,PH_Calibration=Auto,Flow_Order=TACG,Ion_Sphere_Metrics=Good>
##PLATFORM_INDEPENDENT=<Pooling=BatchA,Sequencing_Kit=KitX,Base_Calling_Alg=BaseCaller,Q30=92.5,Normalized_Coverage=120X,Run_Specific_Calibration=Manual>
##BIOINFO_ALIGNMENT=<Software=BWA,Params=-M,Ref_Genome_Version=GRCh38,Recalibration=GATK>
##BIOINFO_VARIANT_CALLING=<Tool=GATK,Version=4.2,Filtering_Thresholds=\"DP>10\",Duplicate_Handling=Marked>
##BIOINFO_POSTPROC=<Normalization=Global,Harmonization=Batch>
##INPUT_QUALITY=<A260_A280=1.8,A260_A230=2.1,DNA_Concentration=15.2,RNA_Concentration=7.1,Notes=\"High quality\",Additional_Metrics=\"{\\\"rna_integrity\\\":7.4,\\\"metrics\\\":{\\\"inner\\\":1}}\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample001
1\t5678\trsMeta\tC\tG\t88\tPASS\tAF=0.25\tGT:GQ\t0/1:80
"""

    MALFORMED_HEADER_VCF = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample001
1\t1234\trsMalformed\tA\tC\t60\tPASS\tSOMATIC;AF=0.1\tGT:GQ\t0/1:45
"""

    FALLBACK_JSON_METADATA_VCF = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=1>\n"
        "##SAMPLE=<ID=JsonGroup,Description=\"JSON fallback\"," 
        "Input_Quality_Additional_Metrics={\"rna_integrity\":7.4,\"metrics\":{\"inner\":1}},"
        "Reference_Genome_Build_Additional_Info={\"notes\":\"fallback\"}>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample001\n"
        "1\t111\trsJson\tA\tG\t50\tPASS\tAF=0.2\tGT:GQ\t0/1:70\n"
    )

    VCF_WITH_COMMON_INFO_FIELDS = """##fileformat=VCFv4.2
##contig=<ID=1>
##INFO=<ID=QD,Number=1,Type=Float,Description="Quality by Depth">
##INFO=<ID=FS,Number=1,Type=Float,Description="Fisher Strand">
##INFO=<ID=SOR,Number=1,Type=Float,Description="Strand Odds Ratio">
##INFO=<ID=CUSTOM,Number=1,Type=String,Description="Custom field">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample001
1\t2000\t.\tA\tG\t.\tPASS\tQD=12.5;FS=7.1;SOR=1.8;CUSTOM=note\tGT:GQ\t0/1:77
"""

    def setUp(self) -> None:
        super().setUp()
        self.user = get_user_model().objects.create_user(
            username="vcf_user",
            email="vcf@example.com",
            password="import-pass",
        )

    def test_import_links_alleles_to_sample_group(self) -> None:
        """Posting a small VCF links the created alleles to the sample group."""

        self.client.login(username="vcf_user", password="import-pass")

        upload = SimpleUploadedFile(
            "import.vcf",
            self.VCF_CONTENT.encode("utf-8"),
            content_type="text/vcf",
        )

        response = self.client.post(reverse("import_data"), {"data_file": upload})

        self.assertRedirects(response, reverse("profile"))

        sample_group = SampleGroup.objects.get(name="GroupA")
        alleles = list(sample_group.allele_frequencies.order_by("id"))

        self.assertEqual(len(alleles), 1)
        allele = alleles[0]

        self.assertEqual(allele.chrom, "1")
        self.assertEqual(allele.pos, 1234)
        self.assertEqual(allele.variant_id, "rsTest")
        self.assertAlmostEqual(float(allele.info.af), 0.5)
        self.assertEqual(allele.info.additional["clnsig"], "Pathogenic")
        self.assertEqual(allele.format.genotype, "0/1")
        self.assertEqual(allele.format.fields["gq"], "99")
        self.assertEqual(allele.format.additional["sample_id"], "Sample001")

        self.assertEqual(
            list(sample_group.allele_frequencies.all()),
            alleles,
        )
        self.assertGreaterEqual(AlleleFrequency.objects.count(), 1)

    def test_import_routes_common_info_fields_to_columns(self) -> None:
        """Common INFO keys map to dedicated columns instead of JSON."""

        self.client.login(username="vcf_user", password="import-pass")

        upload = SimpleUploadedFile(
            "info_fields.vcf",
            self.VCF_WITH_COMMON_INFO_FIELDS.encode("utf-8"),
            content_type="text/vcf",
        )

        response = self.client.post(reverse("import_data"), {"data_file": upload})

        self.assertRedirects(response, reverse("profile"))

        allele = AlleleFrequency.objects.latest("id")

        self.assertAlmostEqual(float(allele.info.qd), 12.5, places=1)
        self.assertAlmostEqual(float(allele.info.fs), 7.1, places=1)
        self.assertAlmostEqual(float(allele.info.sor), 1.8, places=1)
        self.assertEqual(allele.info.additional, {"custom": "note"})

    def test_import_preserves_custom_sample_metadata(self) -> None:
        """Custom SAMPLE header keys are stored and exposed in the detail view."""

        self.client.login(username="vcf_user", password="import-pass")

        upload = SimpleUploadedFile(
            "custom_metadata.vcf",
            self.VCF_CONTENT.encode("utf-8"),
            content_type="text/vcf",
        )

        response = self.client.post(reverse("import_data"), {"data_file": upload})

        self.assertRedirects(response, reverse("profile"))

        sample_group = SampleGroup.objects.get(name="GroupA")
        self.assertEqual(sample_group.additional_metadata, {"customkey": "Value"})

        detail_response = self.client.get(
            reverse("sample_group_detail", args=[sample_group.pk])
        )
        self.assertEqual(detail_response.status_code, 200)
        self.assertContains(detail_response, "Custom metadata")
        self.assertContains(detail_response, "customkey")
        self.assertContains(detail_response, "Value")

    def test_import_populates_related_metadata(self) -> None:
        """Detailed metadata is persisted into structured sample group relations."""

        self.client.login(username="vcf_user", password="import-pass")

        upload = SimpleUploadedFile(
            "metadata.vcf",
            self.VCF_WITH_METADATA.encode("utf-8"),
            content_type="text/vcf",
        )

        response = self.client.post(reverse("import_data"), {"data_file": upload})
        self.assertRedirects(response, reverse("profile"))

        sample_group = SampleGroup.objects.select_related(
            "reference_genome_build",
            "genome_complexity",
            "sample_origin",
            "material_type",
            "bioinfo_alignment",
            "bioinfo_variant_calling",
            "bioinfo_post_proc",
            "input_quality",
        ).get(name="GroupB")

        self.assertEqual(sample_group.comments, "Detailed group")
        self.assertEqual(sample_group.source_lab, "MetaLab")
        self.assertEqual(sample_group.contact_email, "lab@example.com")
        self.assertEqual(sample_group.contact_phone, "+123456789")
        self.assertEqual(sample_group.total_samples, 5)
        self.assertEqual(sample_group.inclusion_criteria, "Adults")
        self.assertEqual(sample_group.exclusion_criteria, "Under18")
        self.assertIsNotNone(sample_group.reference_genome_build)
        self.assertEqual(sample_group.reference_genome_build.build_name, "GRCh38")
        self.assertEqual(sample_group.reference_genome_build.build_version, "v1")
        self.assertEqual(
            sample_group.reference_genome_build.additional_info,
            {"notes": "GRCh"},
        )

        self.assertIsNotNone(sample_group.genome_complexity)
        self.assertEqual(sample_group.genome_complexity.size, "3.2Gb")
        self.assertEqual(sample_group.genome_complexity.ploidy, "Diploid")
        self.assertEqual(sample_group.genome_complexity.gc_content, "41%")

        self.assertIsNotNone(sample_group.sample_origin)
        self.assertEqual(sample_group.sample_origin.tissue, "Blood")
        self.assertEqual(
            sample_group.sample_origin.collection_method,
            "Venipuncture",
        )
        self.assertEqual(
            sample_group.sample_origin.storage_conditions,
            "-80C",
        )
        self.assertEqual(sample_group.sample_origin.time_stored, "6m")

        self.assertIsNotNone(sample_group.material_type)
        self.assertEqual(sample_group.material_type.material_type, "DNA")
        self.assertEqual(sample_group.material_type.integrity_number, "9.5")

        self.assertIsNotNone(sample_group.bioinfo_alignment)
        self.assertEqual(sample_group.bioinfo_alignment.tool, "BWA")
        self.assertEqual(sample_group.bioinfo_alignment.params, "-M")
        self.assertEqual(
            sample_group.bioinfo_alignment.ref_genome_version,
            "GRCh38",
        )
        self.assertEqual(
            sample_group.bioinfo_alignment.recalibration_settings,
            "GATK",
        )

        self.assertIsNotNone(sample_group.bioinfo_variant_calling)
        self.assertEqual(sample_group.bioinfo_variant_calling.tool, "GATK")
        self.assertEqual(sample_group.bioinfo_variant_calling.version, "4.2")
        self.assertEqual(
            sample_group.bioinfo_variant_calling.filtering_thresholds,
            "DP>10",
        )
        self.assertEqual(
            sample_group.bioinfo_variant_calling.duplicate_handling,
            "Marked",
        )

        self.assertIsNotNone(sample_group.bioinfo_post_proc)
        self.assertEqual(sample_group.bioinfo_post_proc.normalization, "Global")
        self.assertEqual(sample_group.bioinfo_post_proc.harmonization, "Batch")

        self.assertIsNotNone(sample_group.input_quality)
        self.assertAlmostEqual(sample_group.input_quality.a260_a280, 1.8)
        self.assertAlmostEqual(sample_group.input_quality.a260_a230, 2.1)
        self.assertAlmostEqual(sample_group.input_quality.dna_concentration, 15.2)
        self.assertAlmostEqual(sample_group.input_quality.rna_concentration, 7.1)
        self.assertEqual(sample_group.input_quality.notes, "High quality")
        self.assertEqual(
            sample_group.input_quality.additional_metrics,
            {"rna_integrity": 7.4, "metrics": {"inner": 1}},
        )

        additional_metadata = sample_group.additional_metadata or {}
        self.assertIn("platform_independent_q30", additional_metadata)
        self.assertAlmostEqual(
            float(additional_metadata["platform_independent_q30"]),
            92.5,
        )

    @mock.patch("app.views.pysam.VariantFile", side_effect=OSError("boom"))
    def test_import_falls_back_to_text_parser(self, mocked_variant_file: mock.Mock) -> None:
        """If pysam fails to open the VCF, the text fallback still imports data."""

        self.client.login(username="vcf_user", password="import-pass")

        upload = SimpleUploadedFile(
            "import.vcf",
            self.VCF_CONTENT.encode("utf-8"),
            content_type="text/vcf",
        )

        response = self.client.post(reverse("import_data"), {"data_file": upload})

        self.assertRedirects(response, reverse("profile"))
        mocked_variant_file.assert_called()

        sample_group = SampleGroup.objects.get(name="GroupA")
        self.assertEqual(sample_group.comments, "Imported group")

        allele = sample_group.allele_frequencies.get()
        self.assertEqual(allele.chrom, "1")
        self.assertEqual(allele.pos, 1234)
        self.assertEqual(allele.variant_id, "rsTest")
        self.assertAlmostEqual(float(allele.info.af), 0.5)
        self.assertEqual(allele.format.genotype, "0/1")
        self.assertEqual(allele.format.additional["sample_id"], "Sample001")

    @mock.patch("app.views.pysam.VariantFile", side_effect=ValueError("malformed header"))
    def test_import_handles_malformed_header_with_text_fallback(
        self, mocked_variant_file: mock.Mock
    ) -> None:
        """A malformed VCF header triggers the text fallback without crashing."""

        self.client.login(username="vcf_user", password="import-pass")

        upload = SimpleUploadedFile(
            "malformed_header.vcf",
            self.MALFORMED_HEADER_VCF.encode("utf-8"),
            content_type="text/vcf",
        )

        response = self.client.post(reverse("import_data"), {"data_file": upload})

        self.assertRedirects(response, reverse("profile"))
        mocked_variant_file.assert_called()

        sample_group = SampleGroup.objects.get()
        # Fallback uses the filename when metadata cannot be extracted from the header
        self.assertEqual(sample_group.name, "malformed_header")

        allele = sample_group.allele_frequencies.get()
        self.assertEqual(allele.chrom, "1")
        self.assertEqual(allele.pos, 1234)
        self.assertEqual(allele.variant_id, "rsMalformed")
        self.assertTrue(allele.info.additional["somatic"])
        self.assertEqual(allele.format.additional["sample_id"], "Sample001")

    @mock.patch("app.views.pysam.VariantFile", side_effect=OSError("fallback"))
    def test_text_fallback_preserves_json_metadata(self, mocked_variant_file: mock.Mock) -> None:
        """JSON-like SAMPLE metadata survives the manual header parsing fallback."""

        self.client.login(username="vcf_user", password="import-pass")

        upload = SimpleUploadedFile(
            "fallback_json.vcf",
            self.FALLBACK_JSON_METADATA_VCF.encode("utf-8"),
            content_type="text/vcf",
        )

        response = self.client.post(reverse("import_data"), {"data_file": upload})

        self.assertRedirects(response, reverse("profile"))
        mocked_variant_file.assert_called()

        sample_group = SampleGroup.objects.select_related(
            "reference_genome_build",
            "input_quality",
        ).get(name="JsonGroup")

        self.assertEqual(sample_group.comments, "JSON fallback")
        self.assertIsNotNone(sample_group.reference_genome_build)
        self.assertEqual(
            sample_group.reference_genome_build.additional_info,
            {"notes": "fallback"},
        )

        self.assertIsNotNone(sample_group.input_quality)
        self.assertEqual(
            sample_group.input_quality.additional_metrics,
            {"rna_integrity": 7.4, "metrics": {"inner": 1}},
        )
