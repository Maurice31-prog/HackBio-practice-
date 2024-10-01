Target identification

A drug target is defined as a biological entity, usually a protein, that can modulate disease phenotypes \[29]. Thus, the identification of prime drug targets is the first and most important step in drug discovery. Conventional drug target identification strategies are performed experimentally, such as identifying differentially expressed genes between normal and diseased cells or tissues and proteins that are highly interconnected with disease-related proteins.

3.1 Experimental approaches

Conventional experimental approaches for target identification require molecular and biochemical studies of disease pathophysiology. Although such studies expand our knowledge of diseases, they can be time-intensive methods for finding promising drug targets. Recently, genome-scale screening technologies, such as haploinsufficiency profiling (HIP), stable isotope labeling by amino acids in cell culture (SILAC), and target deconvolution have been developed to accelerate target identification.

HIP is a genome-wide screening assay for discovering drug targets by sensitizing cells to chemical compounds and identifying gene products associated with the viability of disease cell lines \[30]. The HIP assay is advantageous in that thousands of genes can be evaluated simultaneously and that no prior knowledge of the pathophysiology of the disease is required.

In conventional drug discovery processes, it is often difficult to identify drug targets due to the complicated pathophysiology of many diseases. In this scenario, a reverse strategy can be applied: chemicals capable of modulating disease phenotypes can be screened, and corresponding target proteins can be found \[31]. Diverse methods are employed for target deconvolution, such as protein microarrays, biochemical suppression, and affinity chromatography \[31]. SILAC is an efficient reverse screening strategy that enables the unbiased, comprehensive, and robust identification of target proteins that bind to small molecule probes and drugs \[32,33]. This technique was recently integrated with quantitative mass spectrometry-based proteomics and affinity chromatography, which enables more accurate identification of drug-protein interactions \[32]. Despite the advantages of SILAC, it has several disadvantages that prevent its widespread and practical use: (i) isotope labeling is costly, (ii) sophisticated instruments, such as high-resolution mass spectrometers are required, and (iii) generating chemically immobilized drugs and ensuring their biological activity takes a long time \[34].

3.2 Computational target identification

Experimental approaches are expensive and are generally conducted at low-throughput scale because of their complexity. To overcome these hurdles, in silico methods have been developed to identify potential drug targets \[35]. Target proteins can be computationally predicted from experimental data \[36,37], derived from the literature using text mining \[38], or inferred from protein networks \[39]. Several web servers such as Harmonizome \[40] and the Open Targets Platform \[41] provide lists of potential drug targets predicted using various databases. Alternatively, a reverse docking technique can be used to identify potential protein targets based on the concept that ligands with similar structures may bind to similar proteins with similar binding affinities, displaying similar biological effects \[42–45].

The association-based identification of drug targets is a commonly used approach. For example, the Open Targets Platform \[46] integrates diverse sources, including omics data, experimental results from animal models, and text-mined data from the literature. The platform then ranks genes according to their association with disease \[46]. Several statistical and machine learning-based models, including TarFisDock \[45], TargetHunter \[47], PharmMapper \[48], and Similarity Ensemble Approach \[49], have been developed to predict the biological targets of a queried drug compound (Table 1). Ligand-based protein target discovery is commonly undertaken when no prior knowledge of pathophysiology is available \[50]. Lavecchia reviewed various machine learning models designed to perform ligand searching using molecular descriptors and fingerprints representing the physicochemical properties of a chemical compound \[51–56]. Given that descriptors and fingerprints are a quantitative representation of the chemical and physical characteristics of a compound, they are widely used in the development of predictive models \[57]. A subtractive approach may help refine predicted targets. For example, potential drug targets to treat Helicobacter pylori infection can be identified by removing redundant enzymes, homologous enzymes to those of human or gut flora, extracellular enzymes, non-essential proteins, and other substances from the proteome of H. pylori \[58]. Several freely available protein target databases are listed in Table 2.

3.3 Target validation

Once a target is identified, the next step is to confirm whether the modulation of the biological function of the target affects the disease phenotype \[78]. There are various methods for modulating biological functions and evaluating predicted targets. Of these methods, small interfering RNAs (siRNAs) \[79] are widely used because they can mimic drug effects by repressing translation, resulting in the temporary suppression of the target protein \[80,81]. siRNAs allow investigation of the effects of target inhibition without inhibitors or prior knowledge of the protein structure \[80]. However, for diseases with complex pathophysiology, such as neurological diseases, the degree of repression by siRNAs may affect cellular physiologies differently and may thereby result in contradictory outcomes \[82]. In such cases, animal models in which the target gene is deleted or mutated can be more informative for target confirmation.

View article

Recent advances on computational approach towards potential drug discovery against leishmaniasis

Tushar Joshi, ... Subhash Chandra, in Pathogenesis, Treatment and Prevention of Leishmaniasis, 2021

4.2.1 Target protein

The essential process in the SBDD method is the identification and validation of target proteins (Grant, 2009). The 3D structures of all therapeutically essential proteins are experimentally evaluated by integrative structure biology approaches such as NMR, X-ray crystallography, or cryo-electron microscopy but if a protein structure is not available, the 3D structure of proteins are used to model by in silico methods. Up to date, the total 72 resolve 3D structures of Leishmania spp. have been submitted in RCSB-PDB (Research Collaboratory for Structural Bioinformatics-Protein data bank) database. Many of these proteins can be used as a receptor to find novel drugs against leishmaniasis. Apart from these structures, some essential drug targets of Leishmania species are not present in protein databank. To resolve this problem, homology modeling approach is one of the excellent and most reliable techniques in computational biology as it builds 3D structures of a protein based on the knowledge of 3D structures of homologous proteins having >40% similarity (Song, Lim, & Tong, 2009).

View chapterExplore book

Drug discovery and development

Rohan Palanki, Sourav K. Bose, in Translational Surgery, 2023

Target discovery and target validation

A drug exerts its effects by interacting with proteins, polysaccharides, lipids, and/or nucleic acids in the human body.4 Almost all small molecule drugs and biotherapeutics act by perturbing the function of proteins, as a result of toxicity and low specificity when targeting other macromolecules. Although there are >20,000 human protein-coding genes, not all proteins are suitable for drug interactions and even fewer are appropriate drug targets.10 The druggable proteome is defined as the fraction of proteins that have the ability to bind a drug with high affinity, possess chemical properties complementary to common drug motifs, and are mechanistically linked to disease pathology without significant involvement in critical biological processes.11 If a drug has already been identified for a given target, that target is by definition druggable. However, only ∼500 biomolecules have been successfully targeted to date for therapeutic intervention, dominated by a core set of families (e.g., kinases, ion channels, and G-protein-coupled receptors).12 As such, the development of novel drugs may require the identification of new therapeutic targets and characterization of their underlying molecular mechanisms.

Target discovery is the process by which possible therapeutic targets in a disease process are identified.4 There are several major approaches to target discovery:

•

Genome wide association studies (GWAS) involve scanning large sets of genetic information to statistically associate specific genetic variants, referred to as single nucleotide polymorphisms, with the risk of a disease or its progression. Recently, GWAS analysis was performed to unveil the critical role of the CFH gene (Chromosome 1) in age-related macular degeneration (AMD), the leading cause of blindness in the elderly.13 In 2021, Gemini Therapeutics was granted Food and Drug Administration (FDA) Fast Track designation for GEM103, a recombinant human CFH protein, for its ability to address the clinical progression of AMD in patients carrying CFH loss-of-function variants.

•

Chemical proteomics is a postgenomic version of classical drug affinity chromatography that is coupled to high-resolution mass spectrometry to determine protein targets of existing drugs with unknown or incomplete mechanisms of action.14 The utility of chemical proteomics in target discovery is exemplified by studies of rapamycin, a macrolide immunosuppressant. Despite its known efficacy, the molecular target of this drug was unknown until chemical proteomics revealed its interaction with FKBP12.15 Further investigation of the rapamycin-FKBP12 complex unveiled its inhibition of the mamallian target of rapamycin (mTOR) and subsequently the mTOR pathway.16 Identification of this novel target prompted the development of small molecule derivatives of rapamycin (e.g., everolimus, temsirolimus) with improved pharmacokinetics and reduced immunosuppressive effects that have been approved for the treatment of cancer.17

•

CRISPR-based studies utilize short guide sequences for genetic loci and Cas9 nucleases to induce double-stranded breaks to systematically knockout, inhibit, or activate large numbers of candidate genes. Perturbations that enhance or hinder a disease phenotype can reveal potential drug targets. For example, researchers applied a CRISPR-Cas-based genetic screening approach to a clinical-stage anticancer compound KPT-9274, which had potent activity against various cancer types but no known causal association between its efficacy and cancer cell sensitivity, identifying nicotinamide phosphoribosyltransferase as its molecular target.18

•

In silico predictive models involve identifying potential disease targets through computational analysis of large-scale databases of chemical, pharmacologic, and clinical data. For example, DTINet is a computational pipeline that can predict novel drug-target interactions from a constructed heterogenous network that integrates existing drug-disease-associated information.19 Using this approach, researchers discovered and validated novel interactions between three drugs and three COX proteins.

Target validation is the process of characterizing the function of therapeutic targets.20 In addition to chemical proteomics, there are two major methods by which drug targets are validated:

•

siRNA knockdown involves using small interfering RNA (siRNA) specific to the mRNA transcript of a known genetic target to induce an RNA interference response mediated by the RNA-induced silencing complex. Transient suppression of the target protein product using this method and subsequent observation of the phenotypic effect in vitro can confirm whether drug development against the target is warranted.21

•

Knockout animal models can reproduce the complex interactions that underlie pathophysiological responses to target modulation. To generate a transgenic mouse model, embryonic murine stem cells may be transfected with siRNA expression vectors specific for the desired target, followed by the introduction of cells into tetraploid blastocysts that are implanted into pseudopregnant surrogates. Resultant mice have constitutive knockdown of the protein target, allowing for in vivo functional phenotypic characterization.22

Failure rates during drug development are highest for compounds with a mechanism of action against a previously “undrugged” molecular target and for poorly understood diseases.23 In particular, these drugs suffer from problems related to failures of on-target biological hypotheses and on- and off-target safety concerns. Front-loading of research to target these issues at an early drug discovery stage informs investment of resources into drug development pipelines and reduces the risk of downstream drug failure.

View chapterExplore book

A systematic review of state-of-the-art strategies for machine learning-based protein function prediction

Tian-Ci Yan, ... Tian Xie, in Computers in Biology and Medicine, 2023

1 Introduction

Drug targets are the basis of new drug development and also a difficult point for drug discovery \[1,2]. Therefore, drug target discovery studies are greatly popular in the field of pharmaceutical sciences \[3–5]. Currently, the majority of drug targets are proteins \[6–10]. The classification of protein functions can provide an understanding of the characteristics of similar proteins. Therefore, the functional classification study of drug target proteins can help people to grasp the characteristics of target proteins, which can help to promote the drug target discovery. To conduct functional classification studies of proteins, it is first necessary to annotate protein functions accurately. In the case of the new coronavirus, which is now attacking the world, the study of the function of its host protein will help in the development of effective drugs \[11]. In addition, proteins are the most prominent organic molecules in cells, as they are building blocks of life. One could argue that life cannot exist without proteins. Understanding the functions of proteins is crucial for studying biological systems and their molecular behavior. Thus, precise prediction of protein function is crucial for biotechnology applications, drug development, and human health research \[12–14].

Currently, biochemical investigations such as X-ray diffraction and cryoelectron microscopy generate the most precise and reliable data for protein function annotation \[15,16]. These techniques, however, are typically quite time-consuming, and the supplies and tools used in these tests are very expensive \[17–19]. The number of proteins with the experimentally verified function is substantially lower than the number of newly discovered protein sequences owing to the rapid advancement of genomics technology \[20–22]. Among the 200 million sequences contained in the UniProt database, less than 1% have been confirmed experimentally \[23–28]. Consequently, the need for tests of protein functions cannot be satisfied nowadays by relying solely on experimental techniques. To solve these problems, scientists have begun to predict protein function computationally, which enables high-throughput screening and simultaneous annotation of several proteins, as opposed to experimental validation \[29,30]. In recent years, machine learning has flourished in the field of pharmaceutical research \[31–34] and has gradually become a powerful approach for predicting protein function \[35–38]. To promote improved and more precise protein function predictions, researchers have developed a range of techniques. To improve prediction performance, machine learning algorithms are used, such as convolutional neural network (CNN) \[39–42], k-nearest neighbors (KNN) \[43,44], recurrent neural network (RNN) \[45], and others. Furthermore, integrated algorithms, such as CNN coupled with RNN \[46], three graph neural networks (GNNs) combined \[47], CNN combined with graph convolutional network (GCN) \[48], are also popular.

Based on the various information sources (sequence-based, structure-based, PPI network-based, and multi-information fusion-based) used for protein function prediction, this review is organized into four categories. Fig. 1 shows the flow chart of combining multiple sources of information, such as sequence, structure, and PPI, to finally predict protein function by consensus function. We overview the current state of research into application of machine learning to protein function prediction, noting the existing difficulties and highlighting ground-breaking discoveries, to offer concepts and strategies for future studies.

￼

Sign in to download hi-res image

Fig. 1. Combining structure, sequence and PPI information for protein function prediction. This is an example of combining sequence, structure, and PPI information for protein function prediction. A sequence of unknown function is compared with proteins of known function from the Uniport, STRING, and BioLiP databases to obtain GO predictions, respectively, and finally a consensus function is used to obtain the final prediction.

View article

In silico methods and tools for drug discovery

Bilal Shaker, ... Dokyun Na, in Computers in Biology and Medicine, 2021

3.2 Computational target identification

Experimental approaches are expensive and are generally conducted at low-throughput scale because of their complexity. To overcome these hurdles, in silico methods have been developed to identify potential drug targets \[35]. Target proteins can be computationally predicted from experimental data \[36,37], derived from the literature using text mining \[38], or inferred from protein networks \[39]. Several web servers such as Harmonizome \[40] and the Open Targets Platform \[41] provide lists of potential drug targets predicted using various databases. Alternatively, a reverse docking technique can be used to identify potential protein targets based on the concept that ligands with similar structures may bind to similar proteins with similar binding affinities, displaying similar biological effects \[42–45].

The association-based identification of drug targets is a commonly used approach. For example, the Open Targets Platform \[46] integrates diverse sources, including omics data, experimental results from animal models, and text-mined data from the literature. The platform then ranks genes according to their association with disease \[46]. Several statistical and machine learning-based models, including TarFisDock \[45], TargetHunter \[47], PharmMapper \[48], and Similarity Ensemble Approach \[49], have been developed to predict the biological targets of a queried drug compound (Table 1). Ligand-based protein target discovery is commonly undertaken when no prior knowledge of pathophysiology is available \[50]. Lavecchia reviewed various machine learning models designed to perform ligand searching using molecular descriptors and fingerprints representing the physicochemical properties of a chemical compound \[51–56]. Given that descriptors and fingerprints are a quantitative representation of the chemical and physical characteristics of a compound, they are widely used in the development of predictive models \[57]. A subtractive approach may help refine predicted targets. For example, potential drug targets to treat Helicobacter pylori infection can be identified by removing redundant enzymes, homologous enzymes to those of human or gut flora, extracellular enzymes, non-essential proteins, and other substances from the proteome of H. pylori \[58]. Several freely available protein target databases are listed in Table 2.

View article

Rational Structure-Based Drug Design

Varun Khanna, ... Nikolai Petrovsky, in Encyclopedia of Bioinformatics and Computational Biology, 2019

Rational Structure-Based Drug Design

Choice of the Drug Target and Structure Determination

The choice of the target protein is primarily made based on therapeutic and biological relevance. The ‘druggable’ target is a protein, peptide or nucleic acid whose activity can be modulated by a small molecule (Owens, 2007). Proteins are often used as drug targets due to their significant role in metabolic or signaling pathways specific to a disease condition and can either be activated or inhibited by small molecule drugs. Once the target has been identified it is essential to determine its three-dimensional (3D) structure ideally with a well-defined drug binding pocket to enable high throughput virtual screening in the SBDD approach. The structure of the target can be determined by one of the following methods:

X-ray crystallography and NMR

Both X-ray crystallography and NMR produce data on the relative position of atoms of a molecule. The basis of X-ray crystallography is the scattering of X-rays from electron clouds of atoms whereas NMR measures the interaction of atomic nuclei. The end-product of crystallographic structure determination is an electron density map which is essentially a contour plot indicating positions in the crystal structure where electrons are most likely to be found. This data must be interpreted in terms of a 3D model using semi-automatic computational methodologies. On the other hand, the end-product of an NMR experiment is usually a set of distances between atomic nuclei that define both bonded and non-bonded close contacts in a molecule. These must be interpreted manually to produce a 3D molecular structure using computational tools. The structure determination in each case requires assumptions and approximations, hence the resulting molecular structures obtained may have errors. The choice of technique depends on many factors including molecular weight, ease of solubility and crystallization of the macromolecule under study. X-ray crystallography remains the main workhorse of structure determination for SBDD. Currently, the RCSB Protein Data bank (Parasuraman, 2012) database contains over 130,000 structures of which 90% were solved through X-ray crystallography (Fig. 1).

￼

Sign in to download full-size image

Fig. 1. Incremental increase in the structures of macromolecules in PDB database since year 2000. (Data obtained from PDB).

Homology modeling

In the absence of an experimentally determined 3D structure of a protein, in silico homology modeling can provide structural models that are comparable to the best results achieved experimentally. In general, 30% target-template sequence identity is required to generate a useful structural model (Forrest et al., 2006). This allows researchers to use the generated in silico models for functional analysis and to predict interactions with other molecules. Homology modeling predicts the structure of the target protein primarily by aligning the target sequence called a query with the sequence of one or more known structures called templates and is based on two major assumptions:

1\.

The structure of a protein is encoded in its sequence thus knowing the sequence should at least, in theory, suffice to obtain the structure. This observation was elegantly demonstrated by Anfinsen when he showed that bovine pancreatic ribonuclease, following exposure to a denaturant, could spontaneously regain its native folded structure (Anfinsen, 1973).

2\.

The structure of a protein is more stable and conserved than its sequence. Therefore, closely related similar sequences will essentially adopt similar structures and more distantly related sequences will at least have similar folding, a relationship first identified by Chothia and Lesk (1986).

In practice homology modeling is a multistep process which can be summarized in the following seven steps: a). Template identification, b). Alignment correction, c). Model building d). Loop modeling, e). Side-chain modeling f). Model optimization and g). Model validation. A comprehensive review of homology modeling is beyond the scope of this article. However, interested readers are referred to a series of publications and reviews on homology modeling (Eswar et al., 2006; Fiser and Sali, 2003; Martí-Renom et al., 2000). A list of tools for homology modeling are given in Table 1. Repositories such as Protein model portal (Arnold et al., 2009), Modbase (Pieper et al., 2011) and SWISS-MODEL (Kiefer et al., 2009) contain proteins models generated using various methods.

Table 1. Tools for homology modeling

NameURLSummaryReferenceRaptorXhttp\://raptorx.uchicago.edu/Predict secondary, tertiary, solvent accessibility, disordered regions and binding sites. Remote homology detectionKällberg et al. (2012)Biskithttp\://biskit.pasteur.fr/Python libryary for typical tasks of structural bioinformatics including homology modeling, docking and dynamicsGrünberg et al. (2007)Phyre2http\://www\.sbg.bio.ic.ac.uk/\~phyre2

References 

1\. Zhu W, He X, Cheng K, Zhang L, Chen D, Wang X, et al. Ankylosing spondylitis:

etiology, pathogenesis, and treatments. Bone Res. (2019) 7:22. doi: 10.1038/s41413-019-

0057-8

2\. Alvarez-Navarro C, Lo

́pez de Castro JA. ERAP1 in ankylosing spondylitis:

genetics, biology and pathogenetic role. Curr Opin Rheumatol. (2013) 25:419–25.

doi: 10.1097/BOR.0b013e328362042f

3\. Tsui FW, Haroon N, Reveille JD, Rahman P, Chiu B, Tsui HW, et al. Association

of an ERAP1 ERAP2 haplotype with familial ankylosing spondylitis. Ann Rheum Dis.

(2010) 69:733–6. doi: 10.1136/ard.2008.103804

4\. Evans DM, Spencer CC, Pointon JJ, Su Z, Harvey D, Kochan G, et al. Interaction

between ERAP1 and HLA-B27 in ankylosing spondylitis implicates peptide handling in

the mechanism for HLA-B27 in disease susceptibility. Nat Genet. (2011) 43:761–7.

doi: 10.1038/ng.873

32\. Reveille JD. An update on the contribution of the MHC to AS susceptibility. Clin

Rheumatol. (2014) 33:749–57. doi: 10.1007/s10067-014-2662-7

5\. Mei Y, Pan F, Gao J, Ge R, Duan Z, Zeng Z, et al. Increased serum IL-17 and IL-

23 in the patient with ankylosing spondylitis. Clin Rheumatol. (2011) 30:269–73.

doi: 10.1007/s10067-010-1647-4

6\. Nady S, Ignatz-Hoover J, Shata MT. Interleukin-12 is the optimum cytokine to

expand human Th17 cells in vitro
