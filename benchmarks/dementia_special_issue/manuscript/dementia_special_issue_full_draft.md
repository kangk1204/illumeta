# IlluMeta exposes cohort and neural cell-type dependence in public dementia methylation IDAT reanalysis

Generated draft: 2026-05-08 20:24 Asia/Seoul

## Article type

Full-length research article / methods application manuscript for the Methods special issue "Advancing Mental & Neural Health Assessment - YMETH". The Methods guide for authors lists full-length research articles or reviews as article types; this draft should therefore not be submitted as an application note.

## Abstract
Public Alzheimer disease methylation IDAT datasets can benchmark reproducible EWAS workflows, but secondary reuse often fails at import, quality-control, covariate-governance, cell-composition, and reporting boundaries. We reframed the IlluMeta dementia package as a Methods/YMETH-compatible, cell-type-aware neuroepigenomics workflow evaluation. Six completed AD/dementia public-IDAT cohorts showed strong context dependence: GSE208623, GSE125895, and GSE284764 produced nonzero strict consensus DMPs, whereas GSE134379 and GSE105109 produced no consensus DMPs and GSE66351 produced native-only guarded sensitivity evidence. As a neural cell-type reference axis, GSE306226 Neurons_vs_Microglia produced 6398 strict and 6269 native Minfi/SeSAMe consensus DMPs across 20 neurons and 20 microglia, but its lambda guard triggered. IlluMeta therefore supports auditable public-IDAT reuse by preserving successful signal recovery and negative/guarded evidence, while preventing overclaiming of universal Alzheimer disease methylation biomarkers.

## Keywords

Alzheimer disease; DNA methylation; EWAS; GEO; IDAT; reproducibility; bioinformatics workflow

## Introduction

Dementia epigenome-wide association studies sit at the intersection of biological heterogeneity and workflow heterogeneity. Brain region, cell composition, postmortem interval, age, sex, array generation, and study-specific processing choices can all shape methylation signals before disease effects are tested. Public repositories contain many raw IDAT datasets that could support method benchmarking and secondary neurodegeneration analyses, yet reusing them often requires rebuilding a fragile chain of IDAT import, probe filtering, batch correction, covariate modeling, differential methylation, DMR calling, and reporting.

IlluMeta was developed to make that chain explicit and repeatable. In the present application, we used dementia-related GEO datasets not to claim a new Alzheimer disease biomarker panel, but to test whether the workflow can process raw public IDAT files into interpretable, auditable outputs across different tissues, platforms, and cohort sizes. This distinction is important: a useful methods application paper should show both successful signal recovery and the points where a workflow refuses to overstate noisy or confounded public data.

We therefore chose a staged public-data design. GSE208623 provided a small balanced EPIC blood cohort for rapid proof-of-work. GSE125895 provided a manageable brain validation cohort with multiple regions and covariates. GSE134379 provided a larger brain 450K stress test with 808 arrays. GSE284764 provided a recent EPIC prefrontal-cortex extension with clear AD/control labels. GSE105109 tested a BS-only 5mC-focused branch after excluding oxBS arrays from the same EWAS contrast, and GSE66351 tested an additional frontal/temporal cortex lane while avoiding direct pooling of sorted neuron/glia fractions. The analysis was aligned with the Methods YMETH special issue scope because it evaluates a computational workflow for mental and neural health assessment using public neurodegeneration methylation data, while the interpretation remains conservative because the current study is a public-data reanalysis rather than a new wet-lab validation experiment.

## Methods

### GEO triage and dataset selection

Candidate datasets were identified with IlluMeta GEO search commands using dementia and Alzheimer disease keywords in NCBI GEO [1]. The combined triage table retained IDAT-backed public methylation datasets and recorded GSE identifier, platform, sample count, tissue, diagnosis labels, raw-size estimate, PubMed identifiers where available, and analysis priority. GSE208623 [2,3], GSE125895 [4-6], GSE134379 [7], GSE284764 [8,9], GSE105109 [14,15], and GSE66351 [16,17] were selected for completed analyses because together they cover a rapid peripheral blood proof cohort, multi-region and single-region brain validation lanes, large stress-test cohorts, a recent EPIC prefrontal-cortex extension, a BS-only 5mC-focused brain branch, and an additional frontal/temporal cortex brain branch. GSE197305 was triaged separately because of its size and relevance to the Brains for Dementia Research cohort [10-13], but it was not run because the GEO matrix metadata exposes brain region, sex, and cell-fraction estimates without a sample-level AD/control or neuropathology-stage contrast.

### Full IlluMeta analyses

GSE208623 was analyzed as Alzheimer disease versus control in peripheral blood leukocytes. GSE125895 was analyzed as Alzheimer disease versus control across entorhinal cortex, dorsolateral prefrontal cortex, cerebellum, and hippocampus. GSE134379 was analyzed as Alzheimer disease versus control across cerebellum and middle temporal gyrus. GSE284764 was analyzed as Alzheimer disease versus control in prefrontal cortex; Braak stage was not used as an adjustment variable because it lies on the disease/pathology axis. GSE105109 was analyzed only in the BS assay subset, and GSE66351 was analyzed as a bulk frontal/temporal cortex subset using a manual limma batch-method override after source-level recovery of the large-cohort batch-screening failure. For all completed cohorts, final model designs were audited from branch-specific `*_Metrics.csv` outputs rather than assumed from metadata labels alone. The generated `cross_cohort/adjustment_summary.tsv` records retained metadata covariates, retained and dropped `Cell_Latent*` terms, SVA terms, batch methods, and dropped design terms for each Minfi, strict SeSAMe, and native SeSAMe branch. Cell adjustment is therefore reported as reference-free latent-factor adjustment when RefFreeEWAS was used, not as direct cell-proportion deconvolution unless a reference-based method was available. The GSE125895, GSE134379, GSE284764, GSE105109, and GSE66351 IDAT imports were run under `LC_ALL=C LANG=C` after direct minfi import checks showed that raw IDAT handling was locale-sensitive in this environment.

IlluMeta generated Minfi [18], strict SeSAMe and native SeSAMe [19] branches; branch-specific DMP and DMR outputs; limma-based differential testing [20]; unwanted-variation and batch diagnostics related to SVA concepts [21,22]; DMR outputs grounded in DMRcate-style regional calling [23]; Minfi/SeSAMe consensus DMP tables; branch-comparison metrics; correction adequacy summaries; CRF outputs; methods summaries; and dashboard HTML files. The GSE134379 and GSE66351 runs required additional recovery work because large public cohorts can spend long periods in CPU-bound batch/reporting phases before writing final dashboard artifacts.

### Cross-cohort summarization and figure generation

Cross-cohort tables were generated from completed IlluMeta result directories with `benchmarks/dementia_special_issue/build_full_cross_cohort_summary.py`. The script records sample sizes, significant DMP counts, consensus DMP counts, DMR row counts, branch concordance, pairwise CpG overlap, top consensus CpGs, and top DMRs. Figure 2 was rendered from these tables with `benchmarks/dementia_special_issue/make_full_cross_cohort_figure.py` and exported as PNG, PDF, and SVG with an accompanying source-data workbook, manifest, and DOCX legend. The data-figure validation helper passed for the generated package.

## Results

### GEO triage yielded a staged dementia methylation analysis set

The IlluMeta search identified 24 IDAT-backed dementia-related methylation candidates, and follow-up raw-gate checks separated true IDAT-accessible cohorts from processed-matrix-only candidates. The first completed cohort, GSE208623, contributed 20 Alzheimer disease and 20 control samples from peripheral blood. GSE125895 contributed 82 Alzheimer disease and 187 control samples from four brain regions. GSE134379 contributed 450 Alzheimer disease and 358 control samples from cerebellum and middle temporal gyrus. GSE284764 contributed 119 Alzheimer disease and 98 control samples from prefrontal cortex after mixed-array-size guarding. GSE105109 contributed 136 Alzheimer disease and 56 control samples from the BS-only entorhinal cortex/cerebellum subset. GSE66351 contributed 76 Alzheimer disease and 52 control samples from the frontal/temporal cortex bulk subset. This staged set let us test the workflow on blood, multi-region brain, large brain, EPIC brain, assay-split brain, and additional bulk-cortex contexts.

### Full dual-pipeline outputs were obtained for the proof and validation cohorts

GSE208623 produced 684 hypermethylated and 4,662 hypomethylated Minfi DMPs, 982 and 1,928 strict SeSAMe DMPs, and 1,533 and 2,850 native SeSAMe DMPs. The strict Minfi/SeSAMe consensus contained 2,266 CpGs, and the native consensus contained 3,048 CpGs. DMR calling completed in all three branches, with 4,685, 5,784, and 8,460 rows for Minfi, strict SeSAMe, and native SeSAMe, respectively.

GSE125895 produced a smaller but cleaner brain validation signal. The strict consensus contained 284 CpGs and the native consensus contained 316 CpGs. Branch concordance was high: strict consensus logFC correlation was 0.981 with a significant-set Jaccard overlap of 0.830, and native consensus logFC correlation was 0.950 with a Jaccard overlap of 0.823. DMR outputs again completed in all three branches.

GSE284764 added a recent EPIC prefrontal-cortex brain extension. The strict consensus contained 12 CpGs and the native consensus contained 12 CpGs, with 6,770, 16,168, and 18,327 DMR rows across Minfi, strict SeSAMe, and native SeSAMe. Branch logFC correlation remained high for the nonzero consensus set, but the lambda guard triggered, so this cohort is best framed as guarded supporting evidence rather than a stand-alone positive replication.

GSE105109 and GSE66351 extended the brain 450K lane after raw-IDAT gate failures in other candidates. GSE105109 completed with 1,560, 1,907, and 1,965 DMR rows across Minfi, strict SeSAMe, and native SeSAMe, but produced zero strict and zero native consensus DMPs. GSE66351 completed with 975, 937, and 1,628 DMR rows, zero strict consensus DMPs, and 10 native consensus DMPs. The native GSE66351 branch triggered lambda-guard warning, so those 10 CpGs should be treated as sensitivity evidence rather than a new definitive replication set.

### Cross-cohort consensus overlap was conservative but directionally consistent

The strict and native consensus sets from GSE208623 and GSE125895 shared four CpGs: cg00709979;cg04315214;cg14622549;cg20240347. The overlap is numerically small, which is expected given the tissue and platform differences between peripheral blood EPIC arrays and multi-region brain 450K arrays. The four shared CpGs nevertheless had positive mean logFC in both cohorts, providing a narrow but directionally consistent bridge between the proof-of-work and brain validation analyses.

The brain-to-brain comparison became stronger after adding GSE284764. GSE125895 and GSE284764 shared six strict and native consensus CpGs: cg05066959;cg05810363;cg07584855;cg12309456;cg17104258;cg18102633. These CpGs were directionally concordant across the two brain cohorts. GSE105109 and GSE66351 did not add additional cross-cohort CpG-identity overlaps because GSE105109 had no consensus DMPs and GSE66351's signal was native-only. This does not create a broad Alzheimer disease methylation signature by itself, but it does provide a relevant cross-platform brain replication point while showing where stricter cohort compatibility breaks down. The strongest manuscript claim should therefore emphasize workflow reproducibility, conservative consensus behavior, and small guarded brain consensus bridges rather than a broad disease-signature claim.

### Large, assay-split, and additional brain extensions constrained the strength of the replication claim

GSE134379 generated full branch-level DMP tables and DMR tables, including 3,349 Minfi DMR rows, 3,178 strict SeSAMe DMR rows, and 3,229 native SeSAMe DMR rows. However, it produced only 1 significant Minfi DMP and no strict or native SeSAMe significant DMPs, yielding zero strict and zero native consensus DMPs. The lambda guard triggered. This result is useful because it demonstrates that IlluMeta can finish a large public brain cohort and surface DMR/reference outputs, but it should be handled as a cautionary cohort in the manuscript.

GSE284764 moved the extension back to EPIC arrays and prefrontal cortex. It produced a small but nonzero consensus set, led by CpGs annotated to loci including ANK1/MIR486, RHBDF2, PLVAP, ATG16L2, HOXA3, CTSF, HLX, and DUSP27. Because the lambda guard also triggered, it should be interpreted as guarded support for the workflow's ability to complete a recent raw-IDAT brain EPIC cohort, not as a broad new disease-signature claim. GSE105109 and GSE66351 further narrowed the interpretation: the former was a completed negative consensus cohort with DMR outputs, and the latter showed native-only CpG signal with a native-branch inflation warning. Treating any single extension as a clean positive validation would be incorrect; treating them as stress tests with explicit guardrails is the defensible position.

## Discussion

These analyses show that IlluMeta can move public dementia methylation IDATs through a complete, inspectable analysis path across different cohort sizes and biological contexts. The strongest evidence comes from the combination of a small balanced blood proof cohort, a covariate-aware brain validation cohort, a recent EPIC prefrontal-cortex extension, and two additional brain 450K stress tests that made the boundary conditions visible. The cross-cohort overlap was intentionally evaluated with a conservative CpG-identity criterion: four CpGs were shared between the blood and first brain validation consensus sets, and six CpGs were shared between the 450K brain validation cohort and the EPIC prefrontal-cortex extension. Those overlaps should not be framed as a broad disease signature. They are better interpreted as reproducibility stress tests showing that the workflow can preserve branch agreement while making limited claims when tissue, platform, assay design, and cohort structure differ.

The GSE134379, GSE284764, GSE105109, and GSE66351 extension results are scientifically useful precisely because they constrain the manuscript. A method paper that reports only the successful-looking cohorts would overstate the reliability of public-data reuse. In contrast, GSE134379 shows that a large sample size does not automatically translate into a clean consensus DMP set, GSE284764 shows that a recent EPIC brain cohort can yield a plausible but inflation-guarded consensus set, GSE105109 shows that DMR-level signal can persist without any CpG-level consensus, and GSE66351 shows that native-only signal should not be promoted to a strict consensus claim. IlluMeta's lambda guard and consensus layer provide a practical way to keep those distinctions visible.

Before submission, two constraints should remain visible. First, the Methods guide for authors treats computational work cautiously and asks for experimentally validated biological insight; this manuscript should therefore be framed as an analysis-workflow paper supported by public experimental methylation datasets and independent cohort reanalysis, not as a purely theoretical method. Second, the additional candidate exclusions should be handled transparently. GSE197305 is highly relevant and large, but the sample-level disease or neuropathology contrast was not exposed in GEO metadata. GSE76105, GSE109627, GSE80970, and GSE156984 were raw-blocked for the present raw-IDAT lane, while GSE153712 remains a feasible whole-blood extension that should not be pooled with the brain cohorts. Without external phenotype mapping or new raw-IDAT evidence, the current manuscript should state that six completed cohorts were used and that excluded candidates were blocked by explicit data-access or design criteria rather than ignored.

## Cell-type-aware neuroepigenomics reframe

Generated: 2026-05-22T13:42:38.968915+00:00

### Source-grounded conclusion

The current package is strongest as a Methods/YMETH-compatible neuroepigenomics workflow paper, not as a universal Alzheimer disease biomarker paper. The completed AD/dementia public-IDAT analyses cover 771 controls and 883 AD/test samples across blood, bulk brain, EPIC brain, BS-only brain, and frontal/temporal cortex contexts. They produced 2562 strict and 3386 native consensus DMPs in aggregate, but the evidence is heterogeneous: 3/6 cohorts triggered primary lambda guard; six AD/dementia cohorts are workflow evidence with heterogeneity: 3 strict-zero and 2 native-zero consensus cohorts, plus only 4 nonzero pairwise CpG-overlap rows.

### Guarded AD evidence

- `GSE125895` remains the cleanest brain AD anchor because it has nonzero strict/native consensus and primary lambda guard `ok`.
- `GSE284764` provides recent EPIC prefrontal-cortex support, but its primary lambda guard triggered, so it remains guarded support.
- `GSE134379` and `GSE105109` completed branch-level DMP/DMR outputs but produced no consensus DMPs.
- `GSE66351` should be described as a bulk-cortex sensitivity lane: 0 strict and 10 native consensus DMPs; Bulk frontal/temporal cortex AD run supports heterogeneity/guardrail narrative rather than a universal AD CpG panel.

### Neural cell-type reference axis

`GSE306226/Neurons_vs_Microglia_results` provides the local positive-control/reference axis for cell-type-aware framing. It compares 20 neurons with 20 microglia, passed sample QC, and produced 6398 strict plus 6269 native Minfi/SeSAMe consensus DMPs. Branch concordance was strict logFC r=0.791, Jaccard=0.302; native logFC r=0.761, Jaccard=0.304. Because primary=triggered; Minfi=triggered; SeSAMe=triggered; native=triggered, it should be used as a reference axis demonstrating that IlluMeta can recover a large neural cell-type contrast, not as an unqualified discovery claim.


### Guarded GSE306227 extension

`GSE306227/Neurons_vs_Microglia_results` is retained as a guarded independent sorted-brain extension, not as a core Figure 3 axis in the current package. It has summary-level evidence for 324706 strict and 334911 native consensus DMPs across 18 microglia/control-axis and 19 neuron/test-axis samples, with primary result mode `tier3_ineligible` and lambda guard `triggered`. However, branch metrics present=False, cell summary present=False, and the run log records RefFreeEWAS unavailable=True. Therefore it should be cited only as guarded extension evidence until branch-level design and cell-adjustment artifacts are synced.


### Manuscript claim guardrail

Defensible central claim: IlluMeta exposes context dependence in public neurodegeneration methylation data by making raw-IDAT import, branch-specific covariate/cell/SV adjustment, dual-pipeline agreement, consensus DMP/DMR outputs, and lambda/QC guards auditable. The manuscript should explicitly avoid phrases such as "universal AD biomarker" or "validated AD signature" unless tied to a specific cohort/context and guard status.

## Ethics Statement

This study is a secondary analysis of publicly available, de-identified GEO methylation datasets. No new human participants were recruited and no new biospecimens were collected for this work. The original studies and repositories should be consulted for source-study consent and ethics approvals.

## Data and code availability

The raw public datasets are available from GEO under GSE208623 [3], GSE125895 [6], GSE134379 [7], GSE284764 [9], GSE105109 [15], GSE66351 [17], and the metadata-blocked GSE197305 [13]. All derived analysis artifacts are local to the IlluMeta project workspace. The dataset triage table is `benchmarks/dementia_special_issue/dataset_triage.tsv`. Full result directories are `projects/GSE208623/AD_vs_Control_full_long_results`, `projects/GSE125895/AD_vs_Control_all_regions_results_lcC`, `projects/GSE134379/AD_vs_Control_all_regions_results_lcC`, `projects/GSE284764/AD_vs_Ctl_PFC_results_lcC_run2`, `projects/GSE105109/AD_vs_Control_BS_only_results_lcC`, and `projects/GSE66351/AD_vs_CTRL_bulk_results_lcC_run5_limmavp0_skipcompare`; the neural cell-type reference axis is `projects/GSE306226/Neurons_vs_Microglia_results`. GSE306227 is recorded as a guarded extension in `benchmarks/dementia_special_issue/gse306227_guarded_extension.tsv` until branch-level design and cell-adjustment artifacts are synced. Cross-cohort source tables are stored under `benchmarks/dementia_special_issue/cross_cohort`, including `adjustment_summary.tsv` for branch-level covariate, Cell_Latent, SVA, batch, and dropped-term evidence. The cell-type reframe source table is `benchmarks/dementia_special_issue/cell_type_reframe_summary.tsv`. Figure 2, source data, legend, and manifest are stored under `benchmarks/dementia_special_issue/figures/FIG2`; the cell-type-aware Figure 3 package is stored under `benchmarks/dementia_special_issue/figures/FIG3`.

## Funding

Funding information should be completed by the corresponding author before submission. If no funding supported the work, the Methods guide recommends stating that the research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors.

## Declaration of competing interest

The authors should complete the Elsevier declarations tool before submission. If there are no competing interests, the submission should state that the authors have nothing to declare.

## Declaration of generative AI and AI-assisted technologies in the manuscript preparation process

During preparation of this work, OpenAI Codex was used to support code execution, artifact checking, reference-audit organization, and manuscript drafting. The authors reviewed and edited the generated material and take full responsibility for the content of the manuscript.

## Figure legends

**Figure 2. Cross-cohort IlluMeta evidence from completed public dementia methylation IDAT analyses.** A, sample composition for the six completed cohorts. B, strict and native Minfi/SeSAMe consensus DMP counts, split by direction. C, branch concordance measured by logFC correlation and significant-set Jaccard overlap; zero-consensus cohorts are annotated as no DMP. D, DMR row counts from Minfi, strict SeSAMe, and native SeSAMe branches. E, nonzero pairwise consensus-CpG overlaps across cohorts; zero-overlap pairs remain in the source workbook. F, strict-consensus CpGs shared between the 450K brain validation cohort and the EPIC prefrontal-cortex extension, with mean logFC values from each cohort.

**Figure 3. Cell-type-aware interpretation axis for the IlluMeta dementia reframe.** Consensus DMP counts are shown for the aggregate six-cohort AD/dementia public-IDAT package, the guarded GSE66351 bulk-cortex sensitivity lane, and the GSE306226 neuron-versus-microglia neural cell-type reference axis. The panel shows scale and context dependence and should not be interpreted as a universal Alzheimer disease biomarker figure.

## References

[1] T. Barrett, S.E. Wilhite, P. Ledoux, C. Evangelista, I.F. Kim, M. Tomashevsky, K.A. Marshall, K.H. Phillippy, P.M. Sherman, M. Holko, A. Yefanov, H. Lee, N. Zhang, C.L. Robertson, N. Serova, S. Davis, A. Soboleva, NCBI GEO: archive for functional genomics data sets--update, Nucleic Acids Res. 41 (2013) D991-D995. https://doi.org/10.1093/nar/gks1193.

[2] S. Wu, F. Yang, S. Chao, B. Wang, W. Wang, H. Li, L. Yu, L. He, X. Li, L. Sun, S. Qin, Altered DNA methylome profiles of blood leukocytes in Chinese patients with mild cognitive impairment and Alzheimer's disease, Front. Genet. 14 (2023) 1175864. https://doi.org/10.3389/fgene.2023.1175864.

[3] NCBI Gene Expression Omnibus, Characterization of Global DNA Methylome in Mild Cognitive Impairment and Alzheimer's Disease [dataset], GEO Series GSE208623, 2023. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208623 (accessed 7 May 2026).

[4] S.A. Semick, R.A. Bharadwaj, L. Collado-Torres, R. Tao, J.H. Shin, A. Deep-Soboslay, J.R. Weiss, D.R. Weinberger, T.M. Hyde, J.E. Kleinman, A.E. Jaffe, V.S. Mattay, Integrated DNA methylation and gene expression profiling across multiple brain regions implicate novel genes in Alzheimer's disease, Acta Neuropathol. 137 (2019) 557-569. https://doi.org/10.1007/s00401-019-01966-5.

[5] Z. Li, W. Guo, T. Zeng, J. Yin, K. Feng, T. Huang, Y.D. Cai, Detecting brain structure-specific methylation signatures and rules for Alzheimer's disease, Front. Neurosci. 16 (2022) 895181. https://doi.org/10.3389/fnins.2022.895181.

[6] NCBI Gene Expression Omnibus, Integrated DNA methylation and gene expression profiling across multiple brain regions implicate novel genes in Alzheimer's disease [dataset], GEO Series GSE125895, 2019. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125895 (accessed 7 May 2026).

[7] NCBI Gene Expression Omnibus, Illumina 450K Methylation Data of Alzheimer's Disease: Middle Temporal Gyrus and Cerebellum [dataset], GEO Series GSE134379, 2019. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134379 (accessed 7 May 2026).

[8] V.T. Laroche, R. Cavill, M. Kouhsar, J. Müller, R.A. Reijnders, J. Harvey, A.R. Smith, J. Imm, J. Koetsier, L. Weymouth, L. MacBean, G. Pegoraro, L. Eijssen, B. Creese, G. Kenis, B.M. Tijms, D. van den Hove, K. Lunnon, E. Pishva, Epigenomic subtypes of late-onset Alzheimer's disease reveal distinct microglial signatures, Acta Neuropathol. 151 (2026) 20. https://doi.org/10.1007/s00401-026-02990-y.

[9] NCBI Gene Expression Omnibus, Epigenome-Wide Association Study of the Interaction between Alzheimer's disease and Systemic Infection [dataset], GEO Series GSE284764, 2025. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE284764 (accessed 7 May 2026).

[10] L. Pihlstrøm, G. Shireby, H. Geut, S.P. Henriksen, A.J.M. Rozemuller, J.A. Tunold, E. Hannon, P. Francis, A.J. Thomas, S. Love, J. Mill, W.D.J. van de Berg, M. Toft, Epigenome-wide association study of human frontal cortex identifies differential methylation in Lewy body pathology, Nat. Commun. 13 (2022) 4932. https://doi.org/10.1038/s41467-022-32619-z.

[11] G. Shireby, E.L. Dempster, S. Policicchio, R.G. Smith, E. Pishva, B. Chioza, J.P. Davies, J. Burrage, K. Lunnon, D. Seiler Vellame, S. Love, A. Thomas, K. Brookes, K. Morgan, P. Francis, E. Hannon, J. Mill, DNA methylation signatures of Alzheimer's disease neuropathology in the cortex are primarily driven by variation in non-neuronal cell-types, Nat. Commun. 13 (2022) 5620. https://doi.org/10.1038/s41467-022-33394-7.

[12] F. Grodstein, B. Lemos, J. Yang, K. de Paiva Lopes, R.A. Vialle, N. Seyfried, Y. Wang, G. Shireby, E. Hannon, A. Thomas, K. Brookes, J. Mill, P.L. De Jager, D.A. Bennett, Genetic architecture of epigenetic cortical clock age in brain tissue from older individuals: alterations in CD46 and other loci, Epigenetics 19 (2024) 2392050. https://doi.org/10.1080/15592294.2024.2392050.

[13] NCBI Gene Expression Omnibus, Cortex DNA methylation profiles for the Brains for Dementia Research cohort [dataset], GEO Series GSE197305, 2022. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197305 (accessed 7 May 2026).

[14] A.R. Smith, R.G. Smith, E. Pishva, E. Hannon, J.A. Roubroeks, J. Burrage, C. Troakes, S. Al-Sarraj, C. Sloan, J. Mill, K. Lunnon, Parallel profiling of DNA methylation and hydroxymethylation highlights neuropathology-associated epigenetic variation in Alzheimer's disease, Clin. Epigenetics 11 (2019) 52. https://doi.org/10.1186/s13148-019-0636-y.

[15] NCBI Gene Expression Omnibus, Parallel profiling of DNA methylation and hydroxymethylation in Alzheimer's disease [dataset], GEO Series GSE105109, 2019. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105109 (accessed 7 May 2026).

[16] G. Gasparoni, S. Bultmann, P. Lutsik, T.F.J. Kraus, S. Sordon, J. Vlcek, V. Dietinger, M. Steinmaurer, M. Haider, C.B. Mulholland, T. Arzberger, S. Roeber, M. Riemenschneider, H.A. Kretzschmar, A. Giese, H. Leonhardt, J. Walter, R. Schneider, DNA methylation analysis on purified neurons and glia dissects age and Alzheimer's disease-specific changes in the human cortex, Epigenetics Chromatin 11 (2018) 41. https://doi.org/10.1186/s13072-018-0211-3.

[17] NCBI Gene Expression Omnibus, DNA methylation analysis on purified neurons and glia from human cortex [dataset], GEO Series GSE66351, 2018. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66351 (accessed 7 May 2026).

[18] M.J. Aryee, A.E. Jaffe, H. Corrada-Bravo, C. Ladd-Acosta, A.P. Feinberg, K.D. Hansen, R.A. Irizarry, Minfi: a flexible and comprehensive Bioconductor package for the analysis of Infinium DNA methylation microarrays, Bioinformatics 30 (2014) 1363-1369. https://doi.org/10.1093/bioinformatics/btu049.

[19] W. Zhou, T.J. Triche Jr., P.W. Laird, H. Shen, SeSAMe: reducing artifactual detection of DNA methylation by Infinium BeadChips in genomic deletions, Nucleic Acids Res. 46 (2018) e123. https://doi.org/10.1093/nar/gky691.

[20] M.E. Ritchie, B. Phipson, D. Wu, Y. Hu, C.W. Law, W. Shi, G.K. Smyth, limma powers differential expression analyses for RNA-sequencing and microarray studies, Nucleic Acids Res. 43 (2015) e47. https://doi.org/10.1093/nar/gkv007.

[21] J.T. Leek, W.E. Johnson, H.S. Parker, A.E. Jaffe, J.D. Storey, The sva package for removing batch effects and other unwanted variation in high-throughput experiments, Bioinformatics 28 (2012) 882-883. https://doi.org/10.1093/bioinformatics/bts034.

[22] J.T. Leek, J.D. Storey, Capturing heterogeneity in gene expression studies by surrogate variable analysis, PLoS Genet. 3 (2007) 1724-1735. https://doi.org/10.1371/journal.pgen.0030161.

[23] T.J. Peters, M.J. Buckley, Y. Chen, G.K. Smyth, C.C. Goodnow, S.J. Clark, Calling differentially methylated regions from whole genome bisulphite sequencing with DMRcate, Nucleic Acids Res. 49 (2021) e109. https://doi.org/10.1093/nar/gkab637.
