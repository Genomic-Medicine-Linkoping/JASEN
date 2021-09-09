#!/usr/bin/env nextflow

log.info """\
 J A S E N - P I P E L I N E
 ===================================
 chewbbaca_db_download:   ${params.chewbbaca_db_download}
 local_ariba_db_dir:      ${params.local_ariba_db_dir}
 chewbbacadb_url:         ${params.chewbbacadb_url}
 prodigal_file:           ${params.prodigal_file}
 prodigal_filepath:       ${params.prodigal_filepath}
 sample_ID:               ${params.sample_ID}
 species:                 ${params.species}
 genome_name:             ${params.genome_name}
 input:                   ${params.input}
 outdir:                  ${params.outdir}
 reference:               ${params.reference}
 bwa:                     ${params.bwa}
 ariba_mlst:              ${params.ariba_mlst}
 ariba_pubmlst_name:      ${params.ariba_pubmlst_name}
 aMRFinderPlus_organism:  ${params.aMRFinderPlus_organism}
 aMRFinderPlus_org_name:  ${params.aMRFinderPlus_org_name}
 aMRFinderPlusDb:         ${params.aMRFinderPlusDb}
 """

if (!(params.pkm && params.location)) {
  exit 1, "YOU HAVE TO PROVIDE A LOCATION AND PACKAGE MANAGER PROFILE E.g. 'nextflow run main.nf -profile local,singularity'"
}

process bwa_index_reference{
  label 'min_allocation'

  output:
  file "database.rdy" into bwa_indexes

  """
  if [ ! -f "${params.bwa}/${params.genome_name}.fna.sa" ]; then
    mkdir -p ${params.bwa}/
    cp ${params.reference} ${params.bwa}/
    bwa index ${params.bwa}/${params.genome_name}.fna
  fi
  touch database.rdy
  """
}

process cgmlst_db_init{
  label 'max_allocation'

  output:
  file 'database.rdy' into chewie_init
  file 'chewiedb.zip' into chewie_source

  when:
  params.chewbbaca_db_download && !params.run_tax_analysis

  script:
    """
    export PATH=\$PATH:$baseDir/bin/
    if ${params.chewbbaca_db_download} ; then
    bash cgmlst_db_init.sh -c ${task.cpus} -d ${params.chewbbacadb} -u ${params.chewbbacadb_url}
    else
      touch database.rdy
    fi
    """
}

process ariba_db_download{
  label 'max_allocation'
 
  output:
  file 'database.rdy' into ariba_init
  file 'phenotypes.tsv' into resfinder_phenotypes

  when:
  !params.run_tax_analysis

  """
  if ${params.ariba_db_download} ; then
    ariba getref resfinder resfinder
    ariba prepareref --cdhit_max_memory 0 --force -f ./resfinder.fa -m ./resfinder.tsv --threads ${task.cpus} ${params.aribadb}
    mv resfinder.fa ${params.aribadb}
    mv resfinder.tsv ${params.aribadb}
    wget -vc ${params.resfinder_phenotypes}
    mv phenotypes.txt phenotypes.tsv
    cp phenotypes.tsv ${params.aribadb}/phenotypes.tsv
    touch database.rdy
  else
    touch database.rdy
  fi
  """
}

process ariba_prepare_localdb{
  label 'max_allocation'
 
  output:
  file 'database_local.rdy' into ariba_init_local
  file 'database_non_coding.rdy' into ariba_init_nonc

  when:
  !params.run_tax_analysis

  """
  ariba prepareref --cdhit_max_memory 0 --force --verbose --all_coding yes -f ${params.local_ariba_db_dir}/final/coding.fa --threads ${task.cpus} ${params.aribadb_local}
  cp ${params.local_ariba_db_dir}/final/coding.fa ${params.aribadb_local}
  touch database_local.rdy

  ariba prepareref --cdhit_max_memory 0 --force --verbose --all_coding no -f ${params.local_ariba_db_dir}/final/non-coding.fa --threads ${task.cpus} ${params.aribadb_nonc}
  cp ${params.local_ariba_db_dir}/final/non-coding.fa ${params.aribadb_nonc}
  touch database_non_coding.rdy
  """
}

samples = Channel.fromPath("${params.input}/*.{fastq.gz,fsa.gz,fa.gz,fastq,fsa,fa}")

process fastqc_readqc{
  label 'max_allocation'

  publishDir "${params.outdir}/fastqc", mode: 'copy', overwrite: true

  input:
  file lane1dir from samples

  output:
  file "*_fastqc.html" into fastqc_results

  """
  fastqc ${params.input}/${lane1dir} --format fastq --threads ${task.cpus} -o .
  FILENAME=\$(basename -s .fastq.gz ${params.input}/${lane1dir})
  mkdir -p ${params.outdir}/fastqc
  cp \$FILENAME'_fastqc.zip' ${params.outdir}/fastqc
  """
}

forward = Channel.fromPath("${params.input}/*_R1_*.{fastq.gz,fsa.gz,fa.gz,fastq,fsa,fa}")
reverse = Channel.fromPath("${params.input}/*_R2_*.{fastq.gz,fsa.gz,fa.gz,fastq,fsa,fa}")


process lane_concatination{
  label 'min_allocation'

  publishDir "${params.outdir}/concatinated", mode: 'copy', overwrite: true

  input:
  file 'forward_concat.fastq.gz' from forward.collectFile()
  file 'reverse_concat.fastq.gz' from reverse.collectFile()

  output:
  tuple 'forward_concat.fastq.gz', 'reverse_concat.fastq.gz' into lane_concat

  """
  #Concatination is done via process flow
  """
}

process trimmomatic_trimming{
  label 'max_allocation'
  publishDir "${params.outdir}/trimmomatic", mode: 'copy', overwrite: true

  input:
  tuple forward, reverse from lane_concat

  output:
  tuple "trim_front_pair.fastq.gz", "trim_rev_pair.fastq.gz", "trim_unpair.fastq.gz" into (trimmed_sample_1, trimmed_sample_2, trimmed_sample_3, trimmed_sample_4, trimmed_sample_5, trimmed_sample_6, trimmed_sample_7, trimmed_sample_8)
  file 'trim_out.log' into trimmommatic_out

  """
  trimmomatic PE -threads ${task.cpus} -phred33 ${forward} ${reverse} trim_front_pair.fastq.gz trim_front_unpair.fastq.gz trim_rev_pair.fastq.gz trim_rev_unpair.fastq.gz ILLUMINACLIP:${params.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> trim_out.log
  cat trim_front_unpair.fastq.gz trim_rev_unpair.fastq.gz >> trim_unpair.fastq.gz
  """
}

process ariba_resistancefind{
  label 'max_allocation'
  publishDir "${params.outdir}/ariba", mode: 'copy', overwrite: true, pattern: 'motif_report.tsv'

  input:
  tuple forward, reverse, unpaired from trimmed_sample_4
  file(database_initalization) from ariba_init

  output:
  file 'motif_report.tsv' into ariba_output

  when:
  !params.run_tax_analysis

  """
  ariba run --spades_options careful --gene_nt_extend 50 --assembly_cov 500 --verbose --force --threads ${task.cpus} ${params.aribadb} ${forward} ${reverse} outdir
  cp outdir/report.tsv motif_report.tsv
  """
}

process ariba_resistancefind_local{
  label 'max_allocation'
  publishDir "${params.outdir}/ariba", mode: 'copy', overwrite: true, pattern: 'motif_report_local.tsv'

  input:
  tuple forward, reverse, unpaired from trimmed_sample_5
  file(database_initalization) from ariba_init_local

  output:
  file 'motif_report_local.tsv' into ariba_output_local

  when:
  !params.run_tax_analysis

  """
  ariba run --spades_options careful --gene_nt_extend 50 --assembly_cov 500 --verbose --force --threads ${task.cpus} ${params.aribadb_local} ${forward} ${reverse} outdir
  cp outdir/report.tsv motif_report_local.tsv
  """
}

process ariba_resistancefind_nonc{
  label 'max_allocation'
  publishDir "${params.outdir}/ariba", mode: 'copy', overwrite: true, pattern: 'motif_report_nonc.tsv'

  input:
  tuple forward, reverse, unpaired from trimmed_sample_6
  file(database_initalization) from ariba_init_nonc

  output:
  file 'motif_report_nonc.tsv' into ariba_output_nonc

  when:
  !params.run_tax_analysis

  """
  ariba run --spades_options careful --gene_nt_extend 50 --assembly_cov 500 --verbose --force --threads ${task.cpus} ${params.aribadb_nonc} ${forward} ${reverse} outdir
  cp outdir/report.tsv motif_report_nonc.tsv
  """
}

process ariba_stats{
  label 'min_allocation'

  publishDir "${params.outdir}/ariba", mode: 'copy', overwrite: true
  cpus 1

  input:
  file(report_local) from ariba_output_local
  file(report_nonc) from ariba_output_nonc
  file(report) from ariba_output

  output:
  tuple 'summary.csv', 'motif_report.json', 'motif_report_local.json', 'motif_report_nonc.json' into ariba_summary_output_1
  file 'motif_report.json' into ariba_summary_output_2a
  file 'motif_report_local.json' into ariba_summary_output_2b
  file 'motif_report_nonc.json' into ariba_summary_output_2c

  when:
  !params.run_tax_analysis

  """
  # If any of the files has some contents (more than the tsv header); run ariba summary
  if [[ \$(wc -l <${report}) -ge 2 ]] || [[ \$(wc -l <${report_local}) -ge 2 ]] || [[ \$(wc -l <${report_nonc}) -ge 2 ]]; then
    ariba summary --col_filter n --row_filter n summary ${report} ${report_local} ${report_nonc}
    python3 $baseDir/bin/tsv_to_json.py ${report} motif_report.json
    python3 $baseDir/bin/tsv_to_json.py ${report_local} motif_report_local.json
    python3 $baseDir/bin/tsv_to_json.py ${report_nonc} motif_report_nonc.json
  else
    # Otherwise just create empty json files
    echo -e "#ariba_ref_name\tref_name\tgene\tvar_only\tflag\treads\tcluster\tref_len\tref_base_assembled\tpc_ident\tctg\tctg_len\tctg_cov\tknown_var\tvar_type\tvar_seq_type\tknown_var_change\thas_known_var\tref_ctg_change\tref_ctg_effect\tref_start\tref_end\tref_nt\tctg_start\tctg_end\tctg_nt\tsmtls_total_depth\tsmtls_nts\tsmtls_nts_depth\tvar_description\tfree_text" > motif_report.tsv
    cp motif_report.tsv motif_report_local.tsv
    cp motif_report.tsv motif_report_nonc.tsv
    python3 $baseDir/bin/tsv_to_json.py motif_report.tsv motif_report.json
    python3 $baseDir/bin/tsv_to_json.py motif_report_local.tsv motif_report_local.json
    python3 $baseDir/bin/tsv_to_json.py motif_report_nonc.tsv motif_report_nonc.json
    touch summary.csv
  fi
  """
}

// Make kraken2 DB contents visible for NF
taxok2d = Channel.fromPath("${params.krakendb}/taxo.k2d")
taxodir = Channel.fromPath("${params.krakendb}/**")

process kraken2_decontamination{
  label 'max_allocation'

  publishDir "${params.outdir}/kraken2", mode: 'copy', overwrite: true

  input:
  tuple forward, reverse, unpaired from trimmed_sample_3
  path x from taxok2d
  path y from taxodir

  output:
  tuple "kraken_out.tsv", "kraken_report.tsv" into (kraken2_output, kraken2_output2)

  """
  kraken2 --db ${params.krakendb} --threads ${task.cpus} --output kraken_out.tsv --report kraken_report.tsv --paired ${forward} ${reverse}
  """
}

// Make kaiju DB contents visible for NF
nodes = Channel.fromPath("${params.kaiju_db}/**")

process kaiju_taxonomy{
  label 'max_allocation'

  publishDir "${params.outdir}/kaiju", mode: 'copy', overwrite: true

  input:
  tuple forward, reverse, unpaired from trimmed_sample_8
  path x from nodes

  output:
  tuple "kaiju.out", "kaiju_summary.tsv" into (kaiju_output, kaiju_output2)

  """
  kaiju -z ${task.cpus} -t ${params.kaiju_db}/nodes.dmp -f ${params.kaiju_db}/nr_euk/kaiju_db_nr_euk.fmi -i ${forward} -j ${reverse} -o kaiju.out
  kaiju2table -t ${params.kaiju_db}/nodes.dmp -n ${params.kaiju_db}/names.dmp -r species -o kaiju_summary.tsv kaiju.out -v
  """
}

process spades_assembly{
  label 'max_allocation'

  publishDir "${params.outdir}/spades", mode: 'copy', overwrite: true

  input:
  file(reads) from trimmed_sample_1

  output:
  file 'scaffolds.fasta' into (assembled_sample_1, assembled_sample_2, assembled_sample_3, assembled_sample_4, assembled_sample_5)

  when:
  !params.run_tax_analysis

  script:
  """
  spades.py --threads ${task.cpus} --careful -o . -1 ${reads[0]} -2 ${reads[1]} -s ${reads[2]}
  """
}

process aMRFinderPlus{
  label 'max_allocation'
  publishDir "${params.outdir}/AMRFinderPlus", mode: 'copy', overwrite: true

  input:
  file(scaffolds) from assembled_sample_5

  output:
  file "aMRFinderPlus.tsv" into aMRFinderPlus
  file "IDedRegions.fa"

  when:
  !params.run_tax_analysis

  script:
  if ( params.aMRFinderPlus_organism )
    """
    mkdir -p ${params.aMRFinderPlusDb}
    amrfinder_update -d ${params.aMRFinderPlusDb}
    amrfinder --nucleotide ${scaffolds} \
    --organism ${params.aMRFinderPlus_org_name} \
    -d ${params.aMRFinderPlusDb}/latest \
    -o aMRFinderPlus.tsv \
    --threads ${task.cpus} \
    --nucleotide_output "IDedRegions.fa" \
    --report_all_equal \
    --coverage_min 0.9 \
    --ident_min 0.9 \
    --plus
    """
  else
    """
    mkdir -p ${params.aMRFinderPlusDb}
    amrfinder_update -d ${params.aMRFinderPlusDb}
    amrfinder --nucleotide ${scaffolds} \
    -d ${params.aMRFinderPlusDb}/latest \
    -o aMRFinderPlus.tsv \
    --threads ${task.cpus} \
    --nucleotide_output "IDedRegions.fa" \
    --report_all_equal \
    --coverage_min 0.9 \
    --ident_min 0.9 \
    --plus
    """
}


// TODO: Now same pubmlst db will be created for each run separately but it should be for each new species only
process mlst_lookup{
  label 'max_allocation'

  publishDir "${params.outdir}/mlst", mode: 'copy', overwrite: true

  input:
  tuple forward, reverse, unpaired from trimmed_sample_7
  file contig from assembled_sample_1

  output:
  file 'mlst.json' into mlst_json_output
  file 'mlst_report.tsv' into ariba_mlst

  when:
  !params.run_tax_analysis

  script:
  if ( params.ariba_mlst )
    """
    ariba pubmlstget --verbose '${params.ariba_pubmlst_name}' ${params.ariba_pubmlst_db}
    # Create one level down from ariba mlst calling output dir
    mkdir -p ${params.outdir}/mlst
    ariba run --threads ${task.cpus} ${params.ariba_pubmlst_db}/ref_db $forward $reverse ${params.outdir}/mlst/ariba
    cp ${params.outdir}/mlst/ariba/mlst_report.tsv mlst_report.tsv
    
    mlst $contig --threads ${task.cpus} --json mlst.json --novel novel_mlst.fasta --minid 99.5 --mincov 95
    """
  else
    """
    mlst $contig --threads ${task.cpus} --json mlst.json --novel novel_mlst.fasta --minid 99.5 --mincov 95
    """
}

process chewbbaca_cgmlst{
  label 'max_allocation'
  publishDir "${params.outdir}/cgmlst", mode: 'copy', overwrite: true

  when:
  params.chewbbaca_db_download

  input:
  file contig from assembled_sample_3
  file 'database.dry' from chewie_init

  output:
  file 'cgmlst_alleles.json' into cgmlst_results_a
  file 'cgmlst_stats.json' into cgmlst_results_b

  when:
  !params.run_tax_analysis

  """
  yes | chewBBACA.py AlleleCall --fr -i \${PWD} -g ${params.chewbbacadb}schema --json --cpu ${task.cpus} -o \${PWD} --ptf ${params.prodigal_filepath}
  mv results_*/* .
  mv results_alleles.json cgmlst_alleles.json
  mv results_statistics.json cgmlst_stats.json
  mkdir -p ${params.chewbbacadb}res/
  cp cgmlst_alleles.json ${params.chewbbacadb}res/cgmlst_alleles.json
  cp cgmlst_stats.json ${params.chewbbacadb}res/cgmlst_stats.json
  """
}


process quast_assembly_qc{
  label 'max_allocation'
  publishDir "${params.outdir}/quast", mode: 'copy', overwrite: true

  input:
  file contig from assembled_sample_2

  output:
  file 'report.tsv' into quast_result, quast_result_2

  when:
  !params.run_tax_analysis

  """
  quast.py $contig -o . -r ${params.reference} -t ${task.cpus}
  mkdir -p ${params.outdir}/quast/icarus_viewers/
  cp icarus_viewers/* ${params.outdir}/quast/icarus_viewers/
  cp icarus.html ${params.outdir}/quast/
  cp report.html ${params.outdir}/quast/
  """
}

process quast_json_conversion{  
  label 'min_allocation'  
  publishDir "${params.outdir}/quast", mode: 'copy', overwrite: true
  cpus 1

  input:
  file(quastreport) from quast_result_2

  output:
  file 'quast_report.json' into quast_result_json

  when:
  !params.run_tax_analysis

  """
  python3 $baseDir/bin/quast_to_json.py $quastreport quast_report.json
  """
}


process bwa_read_mapping{
  label 'max_allocation'
  publishDir "${params.outdir}/bwa", mode: 'copy', overwrite: true

  input:
  file(trimmed) from trimmed_sample_2
  file(database_initalization) from bwa_indexes

  output:
  file 'alignment_sorted.bam' into sorted_sample_1, sorted_sample_2

  when:
  !params.run_tax_analysis

  """
  bwa mem -M -t ${task.cpus} ${params.bwa}/${params.genome_name}.fna ${trimmed[0]} ${trimmed[1]} > alignment.sam
  samtools view --threads ${task.cpus} -b -o alignment.bam -T ${params.reference} alignment.sam
  samtools sort --threads ${task.cpus} -o alignment_sorted.bam alignment.bam
  rm -f alignment.bam alignment.sam
  """
}

process samtools_duplicates_stats{
  label 'min_allocation'

  publishDir "${params.outdir}/samtools", mode: 'copy', overwrite: true

  input:
  file(align_sorted) from sorted_sample_1

  output:
  tuple 'samtools_flagstats.txt', 'samtools_total_reads.txt' into samtools_duplicated_results

  when:
  !params.run_tax_analysis

  """
  samtools flagstat ${align_sorted} &> samtools_flagstats.txt
  samtools view -c ${align_sorted} &> samtools_total_reads.txt
  """
}

process picard_markduplicates{
  label 'min_allocation'

  publishDir "${params.outdir}/picard", mode: 'copy', overwrite: true
  cpus 1

  input:
  file(align_sorted) from sorted_sample_2

  output:
  file 'alignment_sorted_rmdup.bam' into deduplicated_sample, deduplicated_sample_2, deduplicated_sample_3
  file 'picard_duplication_stats.txt' into picard_histogram_output

  when:
  !params.run_tax_analysis

  """
  picard MarkDuplicates I=${align_sorted} O=alignment_sorted_rmdup.bam M=picard_duplication_stats.txt REMOVE_DUPLICATES=true
  """
}

process samtools_calling{
  label 'max_allocation'

  publishDir "${params.outdir}/snpcalling", mode: 'copy', overwrite: true

  input:
  file(align_sorted_rmdup) from deduplicated_sample

  output:
  file 'samtools_calls.bam' into called_sample

  when:
  !params.run_tax_analysis

  """
  samtools view -@ ${task.cpus} -h -q 1 -F 4 -F 256 ${align_sorted_rmdup} | grep -v XA:Z | grep -v SA:Z| samtools view -b - > samtools_calls.bam
  """
}


process vcftools_snpcalling{
  label 'min_allocation'

  publishDir "${params.outdir}/snpcalling", mode: 'copy', overwrite: true

  input:
  file(samhits) from called_sample

  output:
  file 'vcftools.recode.bcf' into snpcalling_output

  when:
  !params.run_tax_analysis

  """
  vcffilter="--minQ 30 --thin 50 --minDP 3 --min-meanDP 20"
  bcffilter="GL[0]<-500 & GL[1]=0 & QR/RO>30 & QA/AO>30 & QUAL>5000 & ODDS>1100 & GQ>140 & DP>100 & MQM>59 & SAP<15 & PAIRED>0.9 & EPP>3"

  freebayes -= --pvar 0.7 -j -J --standard-filters -C 6 --min-coverage 30 --ploidy 1 -f ${params.reference} -b ${samhits} -v freebayes.vcf
  bcftools view freebayes.vcf -o unfiltered_bcftools.bcf.gz -O b --exclude-uncalled --types snps
  bcftools index unfiltered_bcftools.bcf.gz
  bcftools view unfiltered_bcftools.bcf.gz -i \${bcffilter} -o bcftools.bcf.gz -O b
  vcftools --bcf bcftools.bcf.gz \${vcffilter} --remove-filtered-all --recode-INFO-all --recode-bcf --out vcftools
  """
}

process snp_translation{
  publishDir "${params.outdir}/snpcalling", mode: 'copy', overwrite: true

  label 'min_allocation'

  input:
  file bcf_file from snpcalling_output

  output:
  tuple 'vcftools.recode.vcf', 'snp_report.tsv' into snp_translated_output
  file 'snp_report.json' into snp_json_output
  file 'bcftools_stats.txt' into bcftools_stats

  when:
  !params.run_tax_analysis

  script:
  """
  bcftools view ${bcf_file} > vcftools.recode.vcf
  gzip -k vcftools.recode.vcf
  bcftools stats -v vcftools.recode.vcf.gz > bcftools_stats.txt
  gatk VariantsToTable -V vcftools.recode.vcf -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F DP -F I16 -F QS -F MQ0F -GF PL -O snp_report.tsv
  python3 $baseDir/bin/tsv_to_json.py snp_report.tsv snp_report.json
  """
}


process picard_qcstats{
  label 'min_allocation'

  publishDir "${params.outdir}/picard", mode: 'copy', overwrite: true

  input:
  file(alignment_sorted_rmdup) from deduplicated_sample_2

  output:
  tuple 'picard_stats.txt', 'picard_insert_distribution.pdf' into picard_output

  when:
  !params.run_tax_analysis

  """
  picard CollectInsertSizeMetrics I=${alignment_sorted_rmdup} O=picard_stats.txt H=picard_insert_distribution.pdf

  """
}

process samtools_deduplicated_stats{
  label 'min_allocation'

  publishDir "${params.outdir}/samtools", mode: 'copy', overwrite: true

  input:
  file(alignment_sorted_rmdup) from deduplicated_sample_3

  output:
  tuple 'samtools_idxstats.tsv', 'samtools_coverage_distribution.tsv', 'samtools_stats.txt' into samtools_deduplicated_output

  when:
  !params.run_tax_analysis

  """
  samtools index ${alignment_sorted_rmdup}
  samtools idxstats ${alignment_sorted_rmdup} &> samtools_idxstats.tsv
  samtools stats --coverage 1,10000,1 ${alignment_sorted_rmdup} > samtools_stats.txt
  grep ^COV samtools_stats.txt | cut -f 2- &> samtools_coverage_distribution.tsv
  """
}

process multiqc_report_taxonomy{
  label 'min_allocation'

  publishDir "${params.outdir}/multiqc_taxonomy", mode: 'copy', overwrite: true

  //More inputs as tracks are added
  input:
    tuple kraken_output, kraken_report from kraken2_output2
    tuple kaiju_raw, kaiju_summary from kaiju_output2

  output:
    file 'multiqc_report.html' into multiqc_output2

  """
  multiqc ${params.outdir} -f -k json -o \$(pwd)
  """
}

process multiqc_report{
  label 'min_allocation'

  publishDir "${params.outdir}/multiqc", mode: 'copy', overwrite: true

  //More inputs as tracks are added
  input:
    file(quast_report) from quast_result
    file(fastqc_report) from fastqc_results
    tuple snp_vcf, snp_tsv from snp_translated_output
    tuple picard_stats, picard_insert_stats from picard_output
    tuple kraken_output, kraken_report from kraken2_output
    tuple kaiju_raw, kaiju_summary from kaiju_output
    tuple samtools_map, samtools_raw from samtools_duplicated_results
    file(bcftools_stats) from bcftools_stats
    file(trimmomatic_log) from trimmommatic_out
    tuple idxstats, cov_dist, stats from samtools_deduplicated_output

  output:
    file 'multiqc_report.html' into multiqc_output1

  when:
  !params.run_tax_analysis

  """
  multiqc ${params.outdir} -f -k json -x ${params.outdir}/mlst/*_out -o \$(pwd)
  """
}

repeats = Channel.fromPath("${params.spa_db}/sparepeats.fasta")
repeat_order = Channel.fromPath("${params.spa_db}/spatypes.txt")

process spa_gene_extraction {
  label 'min_allocation'
  publishDir "${params.outdir}/spaTyper", mode: 'copy', overwrite: true

	input:
		file(scaffolds) from assembled_sample_4
    file(repeats) from repeats
    file(repeat_order) from repeat_order

	output:
    file "spa_${params.sample_ID}.tsv"

	when:
	params.spa_exist && !params.run_tax_analysis

	"""
  spaTyper \
  -f ${scaffolds} \
  -o ${repeat_order} \
  -r ${repeats} \
  --output spa_${params.sample_ID}.tsv
	"""
}

report_Rmd = Channel.fromPath("${params.report_template_file}")
bibliography = Channel.fromPath("${params.bibliography}")
ref_style = Channel.fromPath("${params.reference_style}")

process build_report{
  label 'min_allocation'
  stageInMode "copy"
  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  file (report) from report_Rmd
  file (mlstjson) from mlst_json_output
  file (motif_report) from ariba_summary_output_2a
  file (motif_report_local) from ariba_summary_output_2b
  file (motif_report_nonc) from ariba_summary_output_2c
  file (quastjson) from quast_result_json
  file (snpreport) from snp_json_output
  file (bibliography) from bibliography
  file (phenotypes) from resfinder_phenotypes
  file (ref_style) from ref_style
  file (aMRFinderPlus) from aMRFinderPlus
  file (ariba_mlst_report) from ariba_mlst

  output:
  file("${html_output}")

  when:
  !params.run_tax_analysis

  script:
  html_output = "${params.sample_ID}.html"

  if ( params.chewbbaca_db_download )
    if ( params.spa_exist )
    """
    cp ${params.chewbbacadb}/res/cgmlst_alleles.json cgmlst_alleles.json
    cp ${params.chewbbacadb}/res/cgmlst_stats.json cgmlst_stats.json
    Rscript -e 'rmarkdown::render(input = "${report}", params = list(sample  = "${params.sample_ID}", cgmlst = TRUE, spa_exist = TRUE, spa_path = "${params.outdir}/spaTyper/spa_${params.sample_ID}.tsv"), output_file = "${html_output}")'
    """
    else
    """
    cp ${params.chewbbacadb}/res/cgmlst_alleles.json cgmlst_alleles.json
    cp ${params.chewbbacadb}/res/cgmlst_stats.json cgmlst_stats.json
    Rscript -e 'rmarkdown::render(input = "${report}", params = list(sample  = "${params.sample_ID}", cgmlst = TRUE, spa_exist = FALSE, spa_path = ""), output_file = "${html_output}")'
    """
  else
    if ( params.spa_exist )
    """
    Rscript -e 'rmarkdown::render(input = "${report}", params = list(sample  = "${params.sample_ID}", cgmlst = FALSE, spa_exist = TRUE, spa_path = "${params.outdir}/spaTyper/spa_${params.sample_ID}.tsv"), output_file = "${html_output}")'
    """
    else
    """
    Rscript -e 'rmarkdown::render(input = "${report}", params = list(sample  = "${params.sample_ID}", cgmlst = FALSE, spa_exist = FALSE, spa_path = ""), output_file = "${html_output}")'
    """
}


/* 
 * completion handler
 */
workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> results/$params.sample_ID/$params.sample_ID"+".html\n" : "Oops .. something went wrong" )

	def msg = """\
		Pipeline execution summary
		---------------------------
		Completed at: ${workflow.complete}
		Duration    : ${workflow.duration}
		Success     : ${workflow.success}
		scriptFile  : ${workflow.scriptFile}
		workDir     : ${workflow.workDir}
		exit status : ${workflow.exitStatus}
		errorMessage: ${workflow.errorMessage}
		errorReport :
		"""
		.stripIndent()
	def error = """\
		${workflow.errorReport}
		"""
		.stripIndent()

	logFile = file("${baseDir}.log.complete")
	logFile.text = msg
	logFile.append(error)

  //   def msg = """\
  //     Pipeline execution summary
  //     ---------------------------
  //     Completed at: ${workflow.complete}
  //     Duration    : ${workflow.duration}
  //     Success     : ${workflow.success}
  //     workDir     : ${workflow.workDir}
  //     exit status : ${workflow.exitStatus}
  //     """
  //     .stripIndent()

  // sendMail(to: 'lauri.mesilaakso@regionostergotland.se', subject: 'My pipeline execution', body: msg)
}