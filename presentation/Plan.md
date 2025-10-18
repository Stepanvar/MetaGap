Вот обновленный раздел "Введение: Актуальность проекта" с интегрированной информацией и источниками:

---

## Введение: Актуальность проекта

### Взрывной рост геномных данных
- К 2025 г. ожидается секвенирование до 2 млрд геномов человека — [увеличение на 4–5 порядков за одно десятилетие](https://www.aridhia.com/blog/is-federated-analysis-the-way-forward-for-genomics/#:~:text=Last%20month%2C%20four%20of%20the,1), требуя масштабных систем хранения и анализа
- Геномика превращается в одну из крупнейших областей **«больших данных»**
- [Объём секвенированных данных удваивается каждые 7 месяцев](https://www.weka.io/resources/solution-brief/accelerate-genomic-discovery-petagene-western-digital-and-wekaio/#:~:text=Based%20on%20current%20sequencing%20rates,to%20reach%20over%2040%20exabytes) — темп роста **превышает закон Мура**
- [К 2025 году прогнозируется накопление 2-40 экзабайт](https://www.weka.io/resources/solution-brief/accelerate-genomic-discovery-petagene-western-digital-and-wekaio/#:~:text=Based%20on%20current%20sequencing%20rates,to%20reach%20over%2040%20exabytes) геномных данных только для человека
- [Геномные данные производятся быстрее, чем могут быть осмысленно проанализированы](https://pmc.ncbi.nlm.nih.gov/articles/PMC9154105/#:~:text=genomic%20data%20are%20now%20being,can%20be%20meaningfully%20analyzed%20and)

### Проблема интерпретации редких вариантов
- **Подавляющее число генетических вариантов — очень редкие**, нередко уникальные для отдельных [популяций](https://www.uu.se/en/press/press-releases/2022/2022-05-16-rare-genetic-variants-not-the-major-contributing-factors-to-common-diseases#:~:text=However%2C%20in%20recent%20years%2C%20novel,causing%20genetic%20effects)
- У отдельных индивидов редкие мутации способны **резко повышать риск [заболевания](https://www.uu.se/en/press/press-releases/2022/2022-05-16-rare-genetic-variants-not-the-major-contributing-factors-to-common-diseases#:~:text=This%20suggests%20that%20the%20major,level%20could%20still%20be%20high)**
- [Недавний экспоненциальный рост человеческой популяции привел к избытку редких генетических вариантов](https://pmc.ncbi.nlm.nih.gov/articles/PMC3586590/#:~:text=Recent%20Explosive%20Human%20Population%20Growth,of%20Rare%20Genetic%20Variants)
- В клинической практике редкие варианты часто лежат в основе наследственных болезней, но их интерпретация затруднена из-за скудности сведений в существующих базах данных

### Необходимость точности: влияние SNP на клинический исход
- **Одиночная нуклеотидная замена (SNP) может определять клинический исход** заболевания
- [SNP rs17632542 в гене PSA влияет на функцию белка, прогрессию опухоли и выживаемость при раке простаты](https://www.nature.com/articles/s41467-024-52472-6#:~:text=we%20show%20that%20the%20rs17632542,preclinical%20xenograft%20models%20of%20tumour)
- [SNP rs10251977 связан с чувствительностью к терапии TKI при раке легких](https://ncri.amegroups.org/article/view/4087/html#:~:text=Patients%20carrying%20this%20SNP%20(heterozygous,are%20more%20sensitive%20to%20TKIs)
- [SNP-SNP взаимодействия могут объяснить различия в риске исходов у пациентов с колоректальным раком](https://pmc.ncbi.nlm.nih.gov/articles/PMC9385108/#:~:text=SNP%20interactions%20may%20explain%20the,outcome%20risk%20among%20colorectal%20cancer)
- [SNP в генах DPYD и UGT являются надежными предикторами токсичности химиотерапии](https://pmc.ncbi.nlm.nih.gov/articles/PMC5982750/#:~:text=several%20SNPs%20of%20dihydropyrimidine%20dehydrogenase,mandatory%20in%20their%20detection%20before)
- **Критическая важность:** точность интерпретации вариантов напрямую влияет на жизнь и здоровье пациентов

### Необходимость комплексного подхода
- **Ключевой тезис:** для корректной интерпретации генетических вариантов требуется **комплекс из множества критических параметров**, охватывающих технические, популяционные и организационные аспекты
- **Ни один из существующих биобанков** или референсных генетических ресурсов не обеспечивает полного набора таких параметров
- Без нового подхода остаются «белые пятна» в понимании значимости многих генетических вариантов

### Фрагментация данных и систематические ошибки
- Геномные сведения [хранятся изолированно](https://www.aridhia.com/blog/is-federated-analysis-the-way-forward-for-genomics/#:~:text=Is%20federated%20analysis%20the%20answer%3F) в разных лабораториях по всему миру
- Россия остаётся [недостаточно представленной](https://ouci.dntb.gov.ua/en/works/4EB5Qpql/#:~:text=studies,us%20to%20characterize%20extensive%20genetic) в глобальных референсных базах

**Последствия фрагментации:**
- Снижение **статистической мощности** исследований при малых выборках
- **Дублирование усилий** — разные группы независимо открывают одни и те же варианты
- **Систематические ошибки в анализе** — отсутствие комплексных данных приводит к накоплению ошибок интерпретации
- **Искажение клинической интерпретации** — вариант, распространенный в России, но редкий в Европе, может ошибочно считаться патологически значимым
- [Исследование RUSeq (7452 экзома)](https://ouci.dntb.gov.ua/en/works/4EB5Qpql/#:~:text=compared%20to%20previous%20studies%20allowed,ru) выявило 51 известный патогенный вариант, **значительно чаще встречающийся в РФ, чем в европейских странах**
- Десятки вариантов с серьезным эффектом присутствуют у здоровых российских доноров, [хотя в ClinVar помечены как патогенные](https://ouci.dntb.gov.ua/en/works/4EB5Qpql/#:~:text=compared%20to%20previous%20studies%20allowed,ru)

### Дилемма защиты конфиденциальности и доступности данных
**Законодательные требования:**
- Генетические данные отнесены к чувствительным персональным данным — [GDPR в ЕС строгие правила хранения и обработки](https://medgorod-clinic.ru/stati/eticheskie-i-pravovye-aspekty-geneticheskikh-issledovaniy/#:~:text=1,%D1%81%D0%B2%D0%B5%D0%B4%D0%B5%D0%BD%D0%B8%D0%B9%2C%20%D0%B2%20%D1%82%D0%BE%D0%BC%20%D1%87%D0%B8%D1%81%D0%BB%D0%B5%20%D0%B3%D0%B5%D0%BD%D0%B5%D1%82%D0%B8%D1%87%D0%B5%D1%81%D0%BA%D0%B8%D1%85)
- В России генетическая информация юридически приравнена к биометрическим персональным данным — [152-ФЗ накладывает дополнительные ограничения](https://medgorod-clinic.ru/stati/eticheskie-i-pravovye-aspekty-geneticheskikh-issledovaniy/#:~:text=%D0%92%20%D0%A0%D0%BE%D1%81%D1%81%D0%B8%D0%B8%20%D0%B2%D0%BE%D0%BF%D1%80%D0%BE%D1%81%D1%8B%20%D1%81%D0%BE%D1%85%D1%80%D0%B0%D0%BD%D0%B5%D0%BD%D0%B8%D1%8F%20%D0%BA%D0%BE%D0%BD%D1%84%D0%B8%D0%B4%D0%B5%D0%BD%D1%86%D0%B8%D0%B0%D0%BB%D1%8C%D0%BD%D0%BE%D1%81%D1%82%D0%B8,%D1%85%D1%80%D0%B0%D0%BD%D0%B5%D0%BD%D0%B8%D1%8F%20%D0%B8%20%D0%BF%D0%B5%D1%80%D0%B5%D0%B4%D0%B0%D1%87%D0%B8%20%D0%BF%D0%BE%D0%BB%D1%83%D1%87%D0%B5%D0%BD%D0%BD%D1%8B%D1%85%20%D1%81%D0%B2%D0%B5%D0%B4%D0%B5%D0%BD%D0%B8%D0%B9)
- [Необходима обязательная анонимизация](https://medgorod-clinic.ru/stati/eticheskie-i-pravovye-aspekty-geneticheskikh-issledovaniy/#:~:text=2,%D1%83%D1%82%D0%B5%D1%87%D0%BA%D0%B8%20%D0%B8%20%D0%BD%D0%B5%D0%BF%D1%80%D0%B0%D0%B2%D0%B8%D0%BB%D1%8C%D0%BD%D0%BE%D0%B3%D0%BE%20%D0%B8%D1%81%D0%BF%D0%BE%D0%BB%D1%8C%D0%B7%D0%BE%D0%B2%D0%B0%D0%BD%D0%B8%D1%8F%20%D1%81%D0%B2%D0%B5%D0%B4%D0%B5%D0%BD%D0%B8%D0%B9) для защиты от утечки и неправильного использования сведений

**Проблема баланса:**
- Необходимо одновременно **соблюдать законодательные нормы** (GDPR, 152-ФЗ) **и обеспечивать исследователям доступ к данным**
- Традиционные подходы к обмену данными не обеспечивают этот баланс
- В РФ **пока отсутствует центральная инфраструктура**, которая обеспечивала бы такой баланс

**Решение через FAIR и федеративный анализ:**
- Использование **агрегированной статистики частот аллелей (AF)** вместо обмена сырыми генотипами
- Показатели частоты аллелей не содержат персональных идентификаторов и представляют обобщенные данные по группе людей
- [Федеративный анализ](https://www.aridhia.com/blog/is-federated-analysis-the-way-forward-for-genomics/#:~:text=While%20some%20sharing%20models%20rely,simultaneously%2C%20whereby%20increasing%20research%20efficiency): данные остаются в месте сбора, передаются только суммарные сведения (частоты аллелей)
- Формат [VCF](https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format#:~:text=VCF%2C%20or%20Variant%20Call%20Format%2C,Alliance%20for%20Genomics%20and%20Health) — стандарт для представления вариантов, совместим с существующими инструментами

**Принципы FAIR (Findable, Accessible, Interoperable, Reusable):**
- **Findable:** данные легко обнаружить через [метаданные и уникальные идентификаторы](https://www.tiledb.com/blog/fair-data-principles-explained#:~:text=Findable%3A%20Data%20should%20be%20easy,with%20rich%2C%20machine%2Dactionable%20metadata)
- **Accessible:** данные [доступны через стандартизированные протоколы](https://www.niaid.nih.gov/research/fair-data-principles#:~:text=Accessibility,in%20a%20secure%20and%20straightforward) с контролем доступа
- **Interoperable:** данные совместимы с другими системами через [стандартизированные форматы и онтологии](https://www.tiledb.com/blog/fair-data-principles-explained#:~:text=Interoperable%3A%20Data%20needs%20to%20be,and%20combined)
- **Reusable:** данные хорошо документированы с [четкими лицензиями и провенансом](https://www.niaid.nih.gov/research/fair-data-principles#:~:text=Reusability,adherence%20to%20field%2D%20or%20community)
- [FAIR поддерживает AI/ML анализ](https://www.tiledb.com/blog/fair-data-principles-explained#:~:text=Supports%20AI%20and%20multi%2Dmodal,essential%20for%20scaling%20AI%20and) и мультимодальные данные
- [FAIR обеспечивает воспроизводимость и прослеживаемость](https://www.tiledb.com/blog/fair-data-principles-explained#:~:text=Ensures%20reproducibility%20and%20traceability) результатов
- [FAIR ускоряет открытие лекарств](https://www.tiledb.com/blog/fair-data-principles-explained#:~:text=Scientists%20at%20the%20United%20Kingdom's,from%20a%20few%20weeks%20to) — пример из Oxford Drug Discovery Institute: сокращение времени оценки генов для болезни Альцгеймера с недель до дней

### Необходимость единой платформы с принципами FAIR
- Объединение разрозненных данных в защищённую инфраструктуру позволит [совместно вычислять частоты аллелей](https://ouci.dntb.gov.ua/en/works/4EB5Qpql/#:~:text=ABSTRACT%20Population%20allele%20frequency%20is,information%20between%20major%20medical%20genetic)
- Устранение дублирования и [повышение статистической мощности](https://www.aridhia.com/blog/is-federated-analysis-the-way-forward-for-genomics/#:~:text=Federated%20analysis%20therefore%20promises%20researchers,better%20insight%20and%20drive%20improvement) исследований
- [Интеграция данных разных лабораторий России](https://ouci.dntb.gov.ua/en/works/4EB5Qpql/#:~:text=attempt%20to%20integrate%20genetic%20information,present%20in%20healthy%20donors%20despite) обеспечит представительство российских данных в глобальном контексте
- **Платформа должна обеспечивать баланс между защитой приватности и доступностью данных, придерживаясь принципов FAIR**

***

**Новые источники добавлены:**
- Weka.io — рост данных удваивается каждые 7 месяцев, 40 экзабайт к 2025
- PMC3586590 — экспоненциальный рост популяции и редкие варианты
- PMC9154105 — данные производятся быстрее, чем анализируются
- Nature s41467-024-52472-6 — SNP rs17632542 влияет на клинический исход рака простаты
- NCRI — SNP rs10251977 связан с чувствительностью к TKI
- PMC9385108 — SNP-SNP взаимодействия и риск исхода
- PMC5982750 — SNP DPYD/UGT как предикторы токсичности
- TileDB, NIAID — принципы FAIR и их преимущества
- Примеры применения FAIR в drug discovery

[1](http://medrxiv.org/lookup/doi/10.1101/2021.06.03.21258293)
[2](https://wwwnc.cdc.gov/eid/article/27/2/20-3767_article.htm)
[3](https://dl.acm.org/doi/10.1145/3448016.3457333)
[4](https://dx.plos.org/10.1371/journal.pone.0090422)
[5](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-38)
[6](https://www.semanticscholar.org/paper/13496ba5b21935f5e284337396840f12c4af48aa)
[7](https://www.semanticscholar.org/paper/3489e60a4865f3f2f990572e7d794fd9cb746051)
[8](http://biorxiv.org/lookup/doi/10.1101/103465)
[9](http://ieeexplore.ieee.org/document/5420280/)
[10](https://www.semanticscholar.org/paper/d088e411971f38483b700f62afe1312cddc07105)
[11](https://pmc.ncbi.nlm.nih.gov/articles/PMC3586590/)
[12](https://pmc.ncbi.nlm.nih.gov/articles/PMC10716689/)
[13](https://pmc.ncbi.nlm.nih.gov/articles/PMC6947637/)
[14](https://pmc.ncbi.nlm.nih.gov/articles/PMC6381353/)
[15](https://arxiv.org/pdf/1803.09632.pdf)
[16](https://pmc.ncbi.nlm.nih.gov/articles/PMC5161661/)
[17](https://pmc.ncbi.nlm.nih.gov/articles/PMC9154105/)
[18](https://pmc.ncbi.nlm.nih.gov/articles/PMC11751491/)
[19](https://dx.plos.org/10.1371/journal.pbio.3001669)
[20](https://pmc.ncbi.nlm.nih.gov/articles/PMC4315300/)
[21](https://pmc.ncbi.nlm.nih.gov/articles/PMC4520474/)
[22](https://onlinelibrary.wiley.com/doi/pdfdirect/10.1111/eva.12659)
[23](https://pmc.ncbi.nlm.nih.gov/articles/PMC7318892/)
[24](https://pmc.ncbi.nlm.nih.gov/articles/PMC3596761/)
[25](https://arxiv.org/pdf/2502.16470.pdf)
[26](https://pmc.ncbi.nlm.nih.gov/articles/PMC5721154/)
[27](https://pmc.ncbi.nlm.nih.gov/articles/PMC4668771/)
[28](https://www.biostars.org/p/9575324/)
[29](https://www.youtube.com/watch?v=_cGssGHWaFI)
[30](https://www.weka.io/resources/solution-brief/accelerate-genomic-discovery-petagene-western-digital-and-wekaio/)
[31](https://pmc.ncbi.nlm.nih.gov/articles/PMC8752133/)
[32](https://genestack.com/assets/pdfs/The%20importance%20of%20metadata%20in%20genomics%20and%20the%20FAIR%20principles%20ebook.pdf)
[33](https://pmc.ncbi.nlm.nih.gov/articles/PMC9385108/)
[34](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2023.1145673/full)
[35](https://www.tiledb.com/blog/fair-data-principles-explained)
[36](https://www.nature.com/articles/s41467-024-52472-6)
[37](http://gregoryzynda.com/ncbi/genome/python/2014/03/31/ncbi-genome.html)
[38](https://www.niaid.nih.gov/research/fair-data-principles)
[39](https://pmc.ncbi.nlm.nih.gov/articles/PMC5982750/)
[40](https://academic.oup.com/dnaresearch/article/23/6/517/2647442)
[41](https://www.nature.com/articles/s41597-022-01265-x)
[42](https://ncri.amegroups.org/article/view/4087/html)
[43](https://frontlinegenomics.com/a-guide-to-the-fair-principles-in-biopharma/)
[44](https://www.sciencedirect.com/science/article/pii/S2589004223024100)
[45](https://www.csescienceeditor.org/article/fair-data-what-it-is/)
[46](https://www.sciencedirect.com/science/article/pii/S2643651524000025)
[47](https://www.labvantage.com/blog/fair-data-principles-research-and-manufacturing/)
[48](https://ranchobiosciences.com/2025/01/unlocking-data-potential-understanding-fair-vs-open-data-in-life-sciences/)