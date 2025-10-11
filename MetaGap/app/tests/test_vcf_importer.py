"""Unit tests covering the VCF importer service."""

from __future__ import annotations

import gzip
import shutil
import tempfile
from pathlib import Path
from unittest import mock

from django.contrib.auth import get_user_model
from django.test import TestCase

from ..models import AlleleFrequency, SampleGroup
from ..services.vcf_importer import VCFImporter


class VCFImporterTests(TestCase):
    """Ensure the importer mirrors the previous view-driven behaviour."""

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

    VCF_WITH_ONT_ONLY_METADATA = """##fileformat=VCFv4.2
##contig=<ID=1>
##SAMPLE=<ID=NanoporeOnly,Description=ONT exclusive dataset>
##SEQUENCING_PLATFORM=<Platform=Oxford Nanopore Technologies,Instrument=PromethION,Flow_Cell=FLO-MIN112,Flow_Cell_Version=v1.0,Pore_Type=R10.4,Bias_Voltage=180mV>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample001
1\t2345\trsOnt\tA\tG\t60\tPASS\tAF=0.4\tGT\t0/1
"""

    VCF_WITH_PLATFORM_INDEPENDENT_ONLY_METADATA = """##fileformat=VCFv4.2
##contig=<ID=1>
##SAMPLE=<ID=PlatformNeutral,Description=Platform agnostic dataset>
##PLATFORM_INDEPENDENT=<Run_Name=NeutralRun,Read_Length=100,Device=GenericSeq,Notes=Stored as text>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample001
1\t3456\trsNeutral\tG\tA\t50\tPASS\tAF=0.1\tGT\t0/1
"""

    MALFORMED_HEADER_VCF = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample001
1\t1234\trsMalformed\tA\tC\t60\tPASS\tSOMATIC;AF=0.1\tGT:GQ\t0/1:45
"""

    FALLBACK_JSON_METADATA_VCF = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=1>\n"
        "##SAMPLE=<ID=JsonGroup,Description=\"JSON fallback\"," \
        "Input_Quality_Additional_Metrics={\"rna_integrity\":7.4,\"metrics\":{\"inner\":1}}," \
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

    VCF_WITH_UNKNOWN_FORMAT_FIELD = """##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AB,Number=1,Type=Float,Description="Allele balance">
##SAMPLE=<ID=FormatGroup,Description=Unknown format>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample001
1\t2222\t.\tC\tA\t50\tPASS\t.\tGT:AB\t0/1:0.65
"""

    VCF_WITH_UNKNOWN_SAMPLE_METADATA = """##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##SAMPLE=<ID=MetadataGroup,Platform=AlienSeq,Description=Unknown sample>
##SEQUENCING_PLATFORM=<Platform=AlienSeq,Instrument=CustomSeq>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample001
1\t3333\trsMetadata\tG\tT\t80\tPASS\t.\tGT\t0/1
"""

    VCF_WITHOUT_EXPLICIT_NAME = """##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample001
2\t4444\trsFallback\tT\tC\t42\tPASS\t.\tGT\t0/1
"""

    VCF_WITH_UNSUPPORTED_PLATFORM = """##fileformat=VCFv4.2
##contig=<ID=1>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##SAMPLE=<ID=PlatformGroup,Description=Unsupported platform>
##SEQUENCING_PLATFORM=<Platform=GalacticSeq,Instrument=Galaxy5000,platform_key=extra>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample001
4\t5555\trsPlatform\tA\tC\t90\tPASS\tDP=25;UNDECLARED=foo\tGT:XY\t0/1:baz
"""

    VCF_WITH_UNKNOWN_METADATA_SECTION = """##fileformat=VCFv4.2
##contig=<ID=1>
##SAMPLE_PLATFORM=<ID=BadGroup,Description=Should fail>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample001
1\t123\t.\tA\tT\t50\tPASS\tDP=10\tGT\t0/1
"""

    VCF_WITH_UNDEFINED_INFO_AND_FORMAT = """##fileformat=VCFv4.2
##contig=<ID=1>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##SAMPLE=<ID=UndefinedGroup,Description=Undefined fields>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample001
3\t7777\tundef\tG\tA\t60\tPASS\tDP=.;MISSING=.;UNDECLARED=3\tGT:XY\t0/1:abc
"""

    def setUp(self) -> None:
        super().setUp()
        self.user = get_user_model().objects.create_user(
            username="vcf_user",
            email="vcf@example.com",
            password="import-pass",
        )
        self.temp_dir = Path(tempfile.mkdtemp(prefix="vcf-import"))
        self.addCleanup(lambda: shutil.rmtree(self.temp_dir, ignore_errors=True))

    def _write_vcf(self, content: str, filename: str = "import.vcf") -> Path:
        path = self.temp_dir / filename
        path.write_text(content)
        return path

    def _import(self, content: str, filename: str = "import.vcf") -> tuple[VCFImporter, SampleGroup]:
        path = self._write_vcf(content, filename)
        return self._import_from_path(path)

    def _import_from_path(self, path: Path) -> tuple[VCFImporter, SampleGroup]:
        importer = VCFImporter(self.user)
        sample_group = importer.import_file(str(path))
        return importer, sample_group

    def _refresh_sequencing_platform_relations(self, sample_group: SampleGroup) -> SampleGroup:
        return SampleGroup.objects.select_related(
            "ont_seq",
            "illumina_seq",
            "pacbio_seq",
            "iontorrent_seq",
        ).get(pk=sample_group.pk)

    def _assert_only_ont_metadata(self, sample_group: SampleGroup) -> None:
        sample_group = self._refresh_sequencing_platform_relations(sample_group)

        self.assertIsNone(sample_group.illumina_seq)
        self.assertIsNone(sample_group.pacbio_seq)
        self.assertIsNone(sample_group.iontorrent_seq)
        self.assertIsNotNone(sample_group.ont_seq)

        assert sample_group.ont_seq is not None
        self.assertEqual(sample_group.ont_seq.instrument, "PromethION")
        self.assertEqual(sample_group.ont_seq.flow_cell, "FLO-MIN112")
        self.assertEqual(sample_group.ont_seq.flow_cell_version, "v1.0")
        self.assertEqual(sample_group.ont_seq.pore_type, "R10.4")
        self.assertEqual(sample_group.ont_seq.bias_voltage, "180mV")

    def _assert_platform_independent_metadata(self, sample_group: SampleGroup) -> None:
        sample_group = self._refresh_sequencing_platform_relations(sample_group)

        self.assertIsNone(sample_group.illumina_seq)
        self.assertIsNone(sample_group.ont_seq)
        self.assertIsNone(sample_group.pacbio_seq)
        self.assertIsNone(sample_group.iontorrent_seq)

        metadata = sample_group.additional_metadata or {}
        self.assertIn("platform_independent_run_name", metadata)
        self.assertEqual(metadata["platform_independent_run_name"], "NeutralRun")
        self.assertEqual(metadata.get("platform_independent_read_length"), 100)
        self.assertEqual(metadata.get("platform_independent_device"), "GenericSeq")
        self.assertEqual(
            metadata.get("platform_independent_notes"),
            "Stored as text",
        )

    def test_import_links_alleles_to_sample_group(self) -> None:
        importer, sample_group = self._import(self.VCF_CONTENT)

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
        self.assertGreaterEqual(AlleleFrequency.objects.count(), 1)
        self.assertEqual(importer.warnings, [])

    def test_import_routes_common_info_fields_to_columns(self) -> None:
        importer, _ = self._import(self.VCF_WITH_COMMON_INFO_FIELDS, filename="info_fields.vcf")

        allele = AlleleFrequency.objects.latest("id")

        self.assertAlmostEqual(float(allele.info.qd), 12.5, places=1)
        self.assertAlmostEqual(float(allele.info.fs), 7.1, places=1)
        self.assertAlmostEqual(float(allele.info.sor), 1.8, places=1)
        self.assertEqual(allele.info.additional, {"custom": "note"})
        self.assertEqual(importer.warnings, [])

    def test_import_preserves_custom_sample_metadata(self) -> None:
        importer, sample_group = self._import(self.VCF_CONTENT, filename="custom_metadata.vcf")

        self.assertEqual(sample_group.additional_metadata, {"customkey": "Value"})
        self.assertEqual(importer.warnings, [])

    def test_import_populates_related_metadata(self) -> None:
        importer, sample_group = self._import(self.VCF_WITH_METADATA, filename="metadata.vcf")

        sample_group = SampleGroup.objects.select_related(
            "reference_genome_build",
            "genome_complexity",
            "sample_origin",
            "material_type",
            "bioinfo_alignment",
            "bioinfo_variant_calling",
            "bioinfo_post_proc",
            "input_quality",
        ).get(pk=sample_group.pk)

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
        self.assertEqual(importer.warnings, [])

    def test_import_parses_ont_only_metadata(self) -> None:
        importer, sample_group = self._import(
            self.VCF_WITH_ONT_ONLY_METADATA,
            filename="ont_only.vcf",
        )

        self._assert_only_ont_metadata(sample_group)
        self.assertEqual(importer.warnings, [])

    @mock.patch("app.services.vcf_importer.pysam.VariantFile", side_effect=OSError("ont fallback"))
    def test_text_fallback_parses_ont_only_metadata(
        self, mocked_variant_file: mock.Mock
    ) -> None:
        importer, sample_group = self._import(
            self.VCF_WITH_ONT_ONLY_METADATA,
            filename="ont_only_fallback.vcf",
        )

        self.assertTrue(mocked_variant_file.called)
        self._assert_only_ont_metadata(sample_group)
        self.assertEqual(len(importer.warnings), 1)
        self.assertIn("Falling back to a text parser", importer.warnings[0])

    def test_import_records_platform_independent_only_metadata(self) -> None:
        importer, sample_group = self._import(
            self.VCF_WITH_PLATFORM_INDEPENDENT_ONLY_METADATA,
            filename="platform_independent_only.vcf",
        )

        self._assert_platform_independent_metadata(sample_group)
        self.assertEqual(importer.warnings, [])

    @mock.patch(
        "app.services.vcf_importer.pysam.VariantFile",
        side_effect=OSError("platform independent fallback"),
    )
    def test_text_fallback_records_platform_independent_only_metadata(
        self, mocked_variant_file: mock.Mock
    ) -> None:
        importer, sample_group = self._import(
            self.VCF_WITH_PLATFORM_INDEPENDENT_ONLY_METADATA,
            filename="platform_independent_only_fallback.vcf",
        )

        self.assertTrue(mocked_variant_file.called)
        self._assert_platform_independent_metadata(sample_group)
        self.assertEqual(len(importer.warnings), 1)
        self.assertIn("Falling back to a text parser", importer.warnings[0])

    def test_import_serializes_unknown_format_fields(self) -> None:
        importer, sample_group = self._import(
            self.VCF_WITH_UNKNOWN_FORMAT_FIELD,
            filename="unknown_format.vcf",
        )

        allele = sample_group.allele_frequencies.get()
        self.assertEqual(allele.format.genotype, "0/1")
        self.assertNotIn("ab", allele.format.fields)
        self.assertAlmostEqual(float(allele.format.additional["ab"]), 0.65, places=2)
        self.assertEqual(importer.warnings, [])

    def test_import_records_unhandled_sample_metadata(self) -> None:
        importer, sample_group = self._import(
            self.VCF_WITH_UNKNOWN_SAMPLE_METADATA,
            filename="unknown_metadata.vcf",
        )

        self.assertEqual(sample_group.name, "MetadataGroup")
        self.assertEqual(sample_group.comments, "Unknown sample")

        additional = sample_group.additional_metadata or {}
        self.assertEqual(additional.get("platform"), "AlienSeq")
        self.assertEqual(
            additional.get("platform_independent_platform"),
            "AlienSeq",
        )
        self.assertEqual(
            additional.get("platform_independent_instrument"),
            "CustomSeq",
        )
        self.assertEqual(len(importer.warnings), 1)
        self.assertIn("Unsupported sequencing platform", importer.warnings[0])

    def test_import_without_sample_name_falls_back_to_filename(self) -> None:
        importer, sample_group = self._import(
            self.VCF_WITHOUT_EXPLICIT_NAME,
            filename="fallback_dataset.vcf",
        )

        self.assertEqual(sample_group.name, "fallback_dataset")
        self.assertIsNone(sample_group.comments)
        self.assertIsNone(sample_group.additional_metadata)
        self.assertEqual(importer.warnings, [])

    def test_import_collects_warning_for_unsupported_platform_section(self) -> None:
        importer, sample_group = self._import(
            self.VCF_WITH_UNSUPPORTED_PLATFORM,
            filename="unsupported_platform.vcf",
        )

        allele = sample_group.allele_frequencies.get()
        self.assertEqual(allele.info.dp, "25")
        self.assertEqual(allele.info.additional.get("undeclared"), "foo")
        self.assertEqual(allele.format.additional.get("xy"), "baz")

        self.assertTrue(
            any(
                "Unsupported sequencing platform" in warning
                for warning in importer.warnings
            ),
        )

    def test_import_collects_unknown_metadata_section(self) -> None:
        importer, sample_group = self._import(
            self.VCF_WITH_UNKNOWN_METADATA_SECTION,
            filename="unknown_section.vcf",
        )

        self.assertEqual(sample_group.name, "unknown_section")
        self.assertIsNotNone(sample_group.additional_metadata)
        assert sample_group.additional_metadata is not None
        platform_metadata = sample_group.additional_metadata.get("sample_platform")
        self.assertIsInstance(platform_metadata, dict)
        assert isinstance(platform_metadata, dict)
        self.assertEqual(platform_metadata.get("ID"), "BadGroup")
        self.assertEqual(platform_metadata.get("Description"), "Should fail")
        self.assertIn(
            "Unsupported metadata section 'SAMPLE_PLATFORM'",
            importer.warnings,
        )

    def test_import_handles_undefined_info_and_format_fields(self) -> None:
        importer, sample_group = self._import(
            self.VCF_WITH_UNDEFINED_INFO_AND_FORMAT,
            filename="undefined_fields.vcf",
        )

        allele = sample_group.allele_frequencies.get()
        self.assertIsNone(allele.info.dp)
        self.assertIsNone(allele.info.additional.get("missing"))
        self.assertEqual(allele.info.additional.get("undeclared"), "3")
        self.assertEqual(allele.format.genotype, "0/1")
        self.assertEqual(allele.format.additional.get("xy"), "abc")
        self.assertEqual(importer.warnings, [])

    @mock.patch("app.services.vcf_importer.pysam.VariantFile", side_effect=OSError("boom"))
    def test_import_falls_back_to_text_parser(self, mocked_variant_file: mock.Mock) -> None:
        importer, sample_group = self._import(self.VCF_CONTENT)

        self.assertTrue(mocked_variant_file.called)
        self.assertEqual(sample_group.comments, "Imported group")

        allele = sample_group.allele_frequencies.get()
        self.assertEqual(allele.chrom, "1")
        self.assertEqual(allele.pos, 1234)
        self.assertEqual(allele.variant_id, "rsTest")
        self.assertAlmostEqual(float(allele.info.af), 0.5)
        self.assertEqual(allele.format.genotype, "0/1")
        self.assertEqual(allele.format.additional["sample_id"], "Sample001")
        self.assertEqual(len(importer.warnings), 1)
        self.assertIn("Falling back to a text parser", importer.warnings[0])

    def test_standard_gzip_file_is_inflated_for_pysam(self) -> None:
        path = self.temp_dir / "standard-compression.vcf.gz"
        path.write_bytes(gzip.compress(self.VCF_CONTENT.encode("utf-8")))

        temp_paths: list[Path] = []

        original_inflate = VCFImporter._inflate_gzip_to_temp
        test_case = self

        def tracking_inflate(importer_self: VCFImporter, source_path: str) -> str:
            temp_path = Path(original_inflate(importer_self, source_path))
            temp_paths.append(temp_path)
            test_case.assertTrue(temp_path.exists())
            return str(temp_path)

        importer = VCFImporter(self.user)
        with mock.patch.object(VCFImporter, "_inflate_gzip_to_temp", tracking_inflate):
            sample_group = importer.import_file(str(path))

        self.assertEqual(len(temp_paths), 1)
        self.assertFalse(temp_paths[0].exists())
        self.assertEqual(importer.warnings, [])

        allele = sample_group.allele_frequencies.get()
        self.assertEqual(allele.chrom, "1")
        self.assertEqual(allele.pos, 1234)
        self.assertEqual(allele.variant_id, "rsTest")
        self.assertAlmostEqual(float(allele.info.af), 0.5)
        self.assertEqual(allele.format.genotype, "0/1")
        self.assertEqual(allele.format.additional["sample_id"], "Sample001")

    @mock.patch("app.services.vcf_importer.pysam.VariantFile", side_effect=OSError("gzip fallback"))
    def test_text_fallback_supports_gzip_encoded_vcf(
        self, mocked_variant_file: mock.Mock
    ) -> None:
        path = self.temp_dir / "compressed.vcf.gz"
        path.write_bytes(gzip.compress(self.VCF_CONTENT.encode("utf-8")))

        importer, sample_group = self._import_from_path(path)

        self.assertTrue(mocked_variant_file.called)

        allele = sample_group.allele_frequencies.get()
        self.assertEqual(allele.chrom, "1")
        self.assertEqual(allele.pos, 1234)
        self.assertEqual(allele.variant_id, "rsTest")
        self.assertAlmostEqual(float(allele.info.af), 0.5)
        self.assertEqual(allele.format.genotype, "0/1")
        self.assertEqual(allele.format.additional["sample_id"], "Sample001")
        self.assertEqual(len(importer.warnings), 1)
        self.assertIn("Falling back to a text parser", importer.warnings[0])

    @mock.patch("app.services.vcf_importer.pysam.VariantFile", side_effect=ValueError("malformed header"))
    def test_import_handles_malformed_header_with_text_fallback(
        self, mocked_variant_file: mock.Mock
    ) -> None:
        importer, sample_group = self._import(self.MALFORMED_HEADER_VCF, filename="malformed_header.vcf")

        self.assertTrue(mocked_variant_file.called)
        self.assertEqual(sample_group.name, "malformed_header")

        allele = sample_group.allele_frequencies.get()
        self.assertEqual(allele.chrom, "1")
        self.assertEqual(allele.pos, 1234)
        self.assertEqual(allele.variant_id, "rsMalformed")
        self.assertTrue(allele.info.additional["somatic"])
        self.assertEqual(allele.format.additional["sample_id"], "Sample001")
        self.assertEqual(len(importer.warnings), 1)
        self.assertIn("Falling back to a text parser", importer.warnings[0])

    @mock.patch("app.services.vcf_importer.pysam.VariantFile", side_effect=OSError("fallback"))
    def test_text_fallback_preserves_json_metadata(self, mocked_variant_file: mock.Mock) -> None:
        importer, sample_group = self._import(self.FALLBACK_JSON_METADATA_VCF, filename="fallback_json.vcf")

        self.assertTrue(mocked_variant_file.called)

        sample_group = SampleGroup.objects.select_related(
            "reference_genome_build",
            "input_quality",
        ).get(pk=sample_group.pk)

        self.assertEqual(sample_group.name, "JsonGroup")
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
        self.assertEqual(len(importer.warnings), 1)
        self.assertIn("Falling back to a text parser", importer.warnings[0])
