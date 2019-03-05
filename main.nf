#!/usr/bin/env nextflow

/*
  This pipeline takes a genome assembly and a FASTA protein file containing
  protein translations of one or more gene families to perform prediction with.

  The more protein predictions provided for each gene family member the better,
  as these are clustered and used to generate protein profiles with which
  prediction is performed. Multiple families can be provided as they are
  clustered using BLASTP (via DIAMOND) and MCL.

  Prediction is then performed with AUGUSTUS using the PPX module. A trained
  species is preferred.

  Software Requirements are as follows:
    clojure
    Augustus (tested with 3.3.1 & 3.3.2)
    exonerate*
    diamond*
    mcl*
    perl
    emboss*

    * indicates that conda will handle the installation
    The plan would be to move most/all of these into a singularity container

*/

// To make this more portable, I took some cues from
// https://github.com/nf-core/rnaseq/blob/master/main.nf

def usage() {
  log.info """
  Usage will be
  nextflow run jguhlin/gf-targeted-prediction \
    --assembly "genome.fa" \
    --family "family.fa" \
    --family_name "my_fam" \
    --aug_path "/location/of/augustus" \
    --species_abbrev "Mysp"

  Genes will be named with the species abbreviation and family name
  """
}

params.family = "$baseDir/spermatogenesis_sequences.fa"
params.family_name = "spermatogenesis"
params.assembly = "/Volumes/archive/deardenlab/guhlin/vvul_final/5_funannotate/final.masked.fasta"
params.augustus_path = "/Volumes/archive/deardenlab/guhlin/software/augustus-3.3.1/"
params.species_abbrev = "Vvulg"

family = file(params.family)
family_name = params.family_name
assembly = file(params.assembly)
augpath = file(params.augustus_path)
spabbrev = params.species_abbrev

// Before we do anything, split the assembly to better utilize paralellization

process splitAssembly {
  cache true
  storeDir "./splitAssembly"
  conda 'bioconda::emboss'
  tag { "$assembly" }

  input: file assembly
  output: file("*.fasta") into split_assembly

  """
  seqretsplit -sequence $assembly -outseq seqoutall
  """
}

split_assembly.into{ datasets_exonerate; datasets_blocksearch }

// Generate exonerate hints for providing to augustus later
// This step takes the longest...
process generateExonerateHints {
  storeDir "./exonerate_results"
  conda 'bioconda::exonerate'
  tag { "$seq" }
  input:
    set ID, file(seq) from datasets_exonerate
  output: set file(seq), file("${ID}.hints") into exonerate_hints

  script:
    """
    exonerate --model protein2genome -Q protein -T dna --showtargetgff --query ${params.family} --target ${seq} -S false > ${ID}.gff3
    /Volumes/archive/deardenlab/guhlin/software/augustus-3.3.1/scripts/exonerate2hints.pl --in=${ID}.gff3 --out=${ID}.hints
    """
}

// Generate augustus PPX profiles
// Preprocess splits into families, based off of blastp + mcl
process preprocessProteinFile {
  tag { "${fam}" }
  conda 'bioconda::diamond bioconda::mcl'
  input:
    file fam from family

  output:
    file 'out.mcl_input.I20' into mcl_ch

  script:
  """
  diamond makedb --in $fam --db db
  diamond blastp --query $fam --db db > self-blast.tsv
  cut -f 1,2,12 self-blast.tsv > mcl_input
  mcl mcl_input --abc -I 2
  """
}

// Once families are split into groups, split them, align them, and generate protein profiles
process generateProteinProfiles {
  tag { "${fam}" }
  input:
    file fam from family
    file mcl from mcl_ch

  output:
    file("out*prfl") into prfl_files

  script:
  """
  perl $baseDir/scripts/split_groups.pl $mcl
  """
}

process genePrediction {
  tag { "${hint.baseName}-${prfl.baseName}" }
  publishDir "./predicted_genes_raw"
  input:
    each prfl from prfl_files.collect()
    set file(seq), file(hint) from exonerate_hints

  output:
    file("${hint.baseName}-${prfl.baseName}.gff3") optional true into gp

  script:
    """
    fastBlockSearch --cutoff=0.8 $seq $prfl > ${hint.baseName}-${prfl.baseName}.output
    perl $baseDir/scripts/prfl_process.pl ${hint.baseName}-${prfl.baseName}.output ${prfl} $baseDir/config $seq ${hint.baseName}-${prfl.baseName}.gff3 ${hint}
    """
}

process renameGenes {
  tag {"Rename genes"}
  publishDir "./predicted_genes", mode: 'copy'
  input: file(predictions) from gp.collect()
  output:
    file ("predicted_genes.gff3")
    file ("predicted_genes.aa")

  """
    cat $predictions > predicted_genes.unprocessed.gff3
    clojure $baseDir/scripts/rename_genes.clj $family_name $spabbrev > predicted_genes.gff3
    perl $augpath/scripts/getAnnoFasta.pl predicted_genes.gff3 --seqfile=$assembly
    transeq -sequence predicted_genes3.codingseq -outseq predicted_genes.aa
  """

}
