#!/usr/bin/env nextflow
/*
 * Authors       :
 *
 *
 *      Scott Hazelhurst
 *      Shaun Aron
 *   	Rob Clucas
 *      Eugene de Beste
 *      Lerato Magosi
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2018
 *
 *
 * Description  : Nextflow pipeline for Wits GWAS.
 *
 */

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths




def helps = [ 'help' : 'help' ]

allowed_params = ["input_dir","input_pat","output","output_dir","data","plink_mem_req","covariates","gemma_num_cores","gemma_mem_req","gemma","linear","logistic","chi2","fisher", "work_dir", "scripts", "max_forks", "high_ld_regions_fname", "sexinfo_available", "cut_het_high", "cut_het_low", "cut_diff_miss", "cut_maf", "cut_mind", "cut_geno", "cut_hwe", "pi_hat", "super_pi_hat", "f_lo_male", "f_hi_female", "case_control", "case_control_col", "phenotype", "pheno_col", "batch", "batch_col", "samplesize", "strandreport", "manifest", "idpat", "accessKey", "access-key", "secretKey", "secret-key", "region", "AMI", "instanceType", "instance-type", "bootStorageSize", "boot-storage-size", "maxInstances", "max-instances", "other_mem_req", "sharedStorageMount", "shared-storage-mount", "max_plink_cores", "pheno","big_time","thin"]
GxE_params=['gxe']
allowed_params+=GxE_params

/*JT : append argume boltlmm, bolt_covariates_type */
/*bolt_use_missing_cov --covarUseMissingIndic : “missing indicator method” (via the --covarUseMissingIndic option), which adds indicator variables demarcating missing status as additional covariates. */
ParamBolt=["bolt_ld_scores_col", "bolt_ld_score_file","boltlmm", "bolt_covariates_type",  "bolt_use_missing_cov", "bolt_num_cores", "bolt_mem_req"]
allowed_params+=ParamBolt
ParamFast=["fastlmm","fastlmm_num_cores", "fastlmm_mem_req", "fastlmm_multi", "fastlmmc_bin"]
allowed_params+=ParamFast


params.each { parm ->
  if (! allowed_params.contains(parm.key)) {
    println "\nUnknown parameter : Check parameter <$parm>\n";
  }
}



def params_help = new LinkedHashMap(helps)


params.queue      = 'batch'
params.work_dir   = "$HOME/h3agwas"
params.input_dir  = "${params.work_dir}/input"
params.output_dir = "${params.work_dir}/output"
params.output_testing = "cleaned"
params.thin       = ""
params.covariates = ""
params.chrom      = ""
outfname = params.output_testing



/* Defines the path where any scripts to be executed can be found.
 */



/* Do permutation testing -- 0 for none, otherwise give number */
params.mperm = 1000

/* Adjust for multiple correcttion */
params.adjust = 0

supported_tests = ["chi2","fisher","model","cmh","linear","logistic","boltlmm", "fastlmm"]


params.chi2     = 0
params.fisher   = 0
params.cmh     =  0
params.model   =  0
params.linear   = 0
params.logistic = 0
params.gemma = 0
params.gemma_mem_req = "6GB"
params.gemma_relopt = 1
params.gemma_lmmopt = 4

/*JT Append initialisation variable*/
params.bolt_covariates_type = ""
params.bolt_ld_score_file= ""
params.bolt_ld_scores_col=""
params.boltlmm = 0
params.bolt_num_cores=8
params.bolt_mem_req="6GB"
params.bolt_use_missing_cov=0
/*fastlmm param*/
params.fastlmm = 0
params.fastlmm_num_cores=8
params.fastlmm_mem_req="15GB"
params.fastlmm_multi = 0 
params.fastlmmc_bin =""
/*gxe param*/
param.GxE=""


params.input_pat  = 'raw-GWA-data'

params.sexinfo_available = "false"


params.plink_mem_req = '750MB' // how much plink needs for this
params.other_process_memory = '750MB' // how much other processed need

max_plink_cores = params.max_plink_cores = 4

plink_mem_req = params.plink_mem_req
other_mem_req = params.other_process_memory

params.help = false


data_ch = Channel.fromPath(params.data)

if (params.help) {
    params.each {
    entry ->
      print "Parameter: <$entry.key>    \t Default: $entry.value"
      if (entry.key == 'help')
          println ""
      else {
        help = params_help.get(entry.key)
        if (help)
          print "\n    $help"
        println ""
      }
  }
  System.exit(-1)
}

if(param.GxE!=''){
for(model in c("boltlmm", "fastlmm")){
if(param[model]==1){
print "GxE and "+model+" not compatible"
System.exit(1)
}
}

//---- Modification of variables for pipeline -------------------------------//


def getConfig = {
  all_files = workflow.configFiles.unique()
  text = ""
  all_files.each { fname ->
      base = fname.baseName
      curr = "\n\n*-subsection{*-protect*-url{$base}}@.@@.@*-footnotesize@.@*-begin{verbatim}"
      file(fname).eachLine { String line ->
	if (line.contains("secretKey")) { line = "secretKey='*******'" }
        if (line.contains("accessKey")) { line = "accessKey='*******'" }
        curr = curr + "@.@"+line 
      }
      curr = curr +"@.@*-end{verbatim}\n"
      text = text+curr
  }
  return text
}



// Checks if the file exists
checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n------\nError in your config\nFile $fn does not exist\n\n---\n")
}




bed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString()
bim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString()
fam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString()


gemma_assoc_ch = Channel.create()
/*JT initatilisation of boltlmm_assoc_ch*/
boltlmm_assoc_ch = Channel.create()
fastlmm_assoc_ch = Channel.create()
rel_ch_fastlmm = Channel.create()

pca_in_ch = Channel.create()
assoc_ch  = Channel.create()
raw_src_ch= Channel.create()

Channel
    .from(file(bed),file(bim),file(fam))
    .buffer(size:3)
    .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
    .set { raw_src_ch }


println "\nTesting data            : ${params.input_pat}\n"
println "Testing for phenotypes  : ${params.pheno}\n"
println "Using covariates        : ${params.covariates}\n\n"

if (params.gemma) println "Doing gemma testing"
if (params.chi2) println "Doing chi2 testing"
if (params.linear) println "Doing linear regression testing"
if (params.logistic) println "Doing logistic regression testing"
if (params.fastlmm == 1) println "Doing mixed model with fastlmm "
if (params.boltlmm == 1) println "Doing mixed model with boltlmm "
println "\n"

if (params.thin)
   thin = "--thin ${params.thin}"
else 
   thin = ""

if (params.chrom) 
   chrom = "--chr ${params.chrom}"
else
   chrom = ""

if (thin+chrom) {
  process thin {
    input: 
      set file(bed), file(bim), file(fam) from raw_src_ch
    output:
      /*JT Append initialisation boltlmm_assoc_ch */
      set file("${out}.bed"), file("${out}.bim"), file("${out}.fam") into  ( pca_in_ch, assoc_ch, gemma_assoc_ch, boltlmm_assoc_ch,fastlmm_assoc_ch, rel_ch_fastlmm)
    script:
       base = bed.baseName
       out  = base+"_t"
       "plink --keep-allele-order --bfile $base $thin $chrom --make-bed --out $out"
  }

  println "\nData has been thinned or only some chromosomes used  (is the run for test purposes only?)\n"
   


} else {
    /*JT : append boltlmm_assoc_ch and a]*/
    raw_src_ch.separate( pca_in_ch, assoc_ch, gemma_assoc_ch, boltlmm_assoc_ch, fastlmm_assoc_ch,rel_ch_fastlmm) { a -> [a,a,a,a,a,a] }
}





process computePCA {
  cpus max_plink_cores
  memory plink_mem_req
  input:
    set file('cleaned.bed'),file('cleaned.bim'),file('cleaned.fam') from pca_in_ch

  publishDir params.output_dir, overwrite:true, mode:'copy'
  output:
    set file("${outfname}.eigenval"), file("${outfname}.eigenvec")  \
         into pca_out_ch
  script:
      base = "cleaned"
      prune= "${base}-prune"
     """
     plink --bfile ${base} --indep-pairwise 100 20 0.2 --out check
     plink --keep-allele-order --bfile ${base} --extract check.prune.in --make-bed --out $prune
     plink --threads $max_plink_cores --bfile $prune --pca --out ${outfname}
     """
}

process drawPCA {
    input:
      set file(eigvals), file(eigvecs) from pca_out_ch
    output:
      file (output) into report_pca_ch
    publishDir params.output_dir, overwrite:true, mode:'copy',pattern: "*.pdf"
    script:
      base=eigvals.baseName
      cc_fname = 0
      cc       = 0
      col      = 0
      // also relies on "col" defined above
      output="${base}-pca.pdf"
      template "drawPCA.py"

}





num_assoc_cores = params.mperm == 0 ? 1 : Math.min(10,params.max_plink_cores)

supported_tests = ["chi2","fisher","model","cmh","linear","logistic"]

requested_tests = supported_tests.findAll { entry -> params.get(entry) }


covariate = ""
gotcovar  = 0
pheno     = ""



 
if (params.data != "") {

   checker(file(params.data))

   if (params.covariates != "") {
      gotcovar = 1
  }

  data_ch1 = Channel.create()
  data_ch2 = Channel.create()
  Channel.fromPath(params.data).separate(data_ch1,data_ch2) { a -> [a,a] } 
  
  process extractPheno {
    input:
     file(data) from data_ch1
    output:
     file(phenof) into pheno_ch
    script:
     phenof = "pheno.phe"
     all_phenos = params.covariates.length()>0 ? params.pheno+","+params.covariates : params.pheno
     """
     extractPheno.py $data ${all_phenos} $phenof
     """
  }


  pheno_label_ch = Channel.from(params.pheno.split(","))

  process showPhenoDistrib {
    input:
    file(data) from data_ch2
    output:
      file ("B050*") into report_ch
    script:
      "phe_distrib.py --pheno ${params.pheno} $data B050 "
  }
}  else {
  report_ch = Channel.empty()
  pheno_label = ""
  pheno_label_ch = Channel.from("")
}

/*JT : Case fastlmm => if yes*/
if (params.fastlmm == 1) {
  data_ch_fastlmm = Channel.fromPath(params.data)
  if(params.fastlmmc_bin=="")fastlmmc="fastlmmc"
  else fastlmmc=params.fastlmmc_bin

  fam_ch_fast = Channel.create()
  gem_ch_fast2 = Channel.create()
  gem_ch_fast =Channel.create()
  bim_ch_fast = Channel.create()
  fastlmm_assoc_ch.separate (gem_ch_fast2,gem_ch_fast,bim_ch_fast,fam_ch_fast) { a -> [a,a, a[1],a[2]] }


  if (params.covariates)
     covariate_option = "--cov_list ${params.covariates}"
  else
     covariate_option = ""

  process  getFastLmmPhenosCovar {
    input:
      file(covariates) from data_ch_fastlmm
      file(fam) from fam_ch_fast
    output:
      set file(phef), file(covfile) into fastlmm_data_ch
      stdout into pheno_cols_ch_fastlmm
    script:
      base = fam.baseName
      phef = "${base}_fastlmm_n.phe"
      covfile = "${base}_fastlmm_n.cov"
      """
      all_covariate.py --data  $covariates --inp_fam  $fam $covariate_option \
                          --pheno ${params.pheno} --phe_out ${phef}  --cov_out $covfile --form_out 3
      """
  }

  ind_pheno_cols_ch = Channel.create()
  check = Channel.create()
  pheno_cols_ch_fastlmm.flatMap { list_str -> list_str.split() }.tap ( check) .set { ind_pheno_cols_ch }

  if(params.fastlmm_multi==1){

     process getRelForFastLMM {
	cpus params.fastlmm_num_cores
	memory params.fastlmm_mem_req
	time params.big_time
	input:
	   file plinks from rel_ch_fastlmm
	output:
	   file("output/${base}.*XX.txt")
	   file("${rel_fastlmm}") into rel_mat_ch_fastlmm
	script:
	   base = plinks[0].baseName
	   fam = plinks[2]
	   rel_fastlmm="rel_fastlmm.txt"
	   """
           cat $fam |awk '{print \$1"\t"\$2"\t"0.2}' > pheno
	   export OPENBLAS_NUM_THREADS=${params.fastlmm_num_cores}
	   gemma -bfile $base  -gk ${params.gemma_relopt} -o $base -p pheno -n 3
	   cvt_rel_gemma_fastlmm.py $fam output/${base}.*XX.txt $rel_fastlmm
	   """
	 }


     process getListeChro{
	input :
	  file(BimFile) from bim_ch_fast
	output :
	  stdout into (chrolist,chrolist2)
	script:
	 """
	 cat $BimFile|awk '{print \$1}'|uniq|sort|uniq
	"""
     }

     check2 = Channel.create()
     ListeChro2=chrolist.flatMap { list_str -> list_str.split() }.tap ( check2)


     process doFastlmmMulti{
       cpus params.fastlmm_num_cores
       time   params.big_time
       input:
	 set file (phef), file(covariate) from fastlmm_data_ch
	 file(rel) from rel_mat_ch_fastlmm
	 file(plinks) from  gem_ch_fast
       each this_pheno from ind_pheno_cols_ch
       each chro from ListeChro2
       output:
	 set (our_pheno, file("$out"), val(base)) into (fastlmm_manhatten_chro,fastlmm_manhatten_chro2)
       script:
	 base = plinks[0].baseName
	 our_pheno = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/-$/,"")
	 covar_opt_fast =  (params.covariates) ?  " -covar $covariate" : ""
	 newbase=base+"-"+chro
	 out = "$base-$our_pheno"+"-"+chro+".stat"
	 """
	 this_pheno_col=`echo ${this_pheno} | sed 's/-.*//' `
	 plink --keep-allele-order --bfile $base --chr $chro --make-bed --out $newbase
	 $fastlmmc -REML -simType RRM -verboseOut -sim $rel -bfile $newbase -pheno ${phef} -simLearnType Full -out $out -maxThreads $params.fastlmm_num_cores \
	          $covar_opt_fast  -mpheno \${this_pheno_col} -bfileSim $base
	 """
       }

     fastlmm_manhatten_chroM=fastlmm_manhatten_chro.groupTuple()
     fastlmm_manhatten_chroM1=fastlmm_manhatten_chro2.groupTuple()


     process doMergeFastlmm{
          input :
	    set (val(this_pheno),list_file, base_list) from fastlmm_manhatten_chroM
	    /* with uniq channels vs 2 => problems*/
	    /*file(plinks) from  gem_ch_fast2*/
	 publishDir "${params.output_dir}/fastlmm", overwrite:true, mode:'copy'
	 output :
	     set val(base), val(our_pheno2), file("$out") into fastlmm_manhatten_ch
	 script :
	     base=base_list[0]
	     our_pheno = this_pheno.replace(/_|\/np.\w+/,"-").replace(/-$/,"")
	     our_pheno2 = this_pheno.replace(/_|\/np.\w+/,"-").replace(/-$/,"").replaceAll(/^[0-9]+-/,"")
	     out = "$base-${our_pheno}.stat"
	     fnames = list_file.join(" ")
	     file1  = list_file[0]
	     """
	     head -1 $file1 > $out
	     cat $fnames | grep -v "Chromosome" >> $out
	     """
     }
  }  else { // if not   doing fastlmm_multi

     process doFastlmm{
       cpus params.fastlmm_num_cores
       time   params.big_time
       input:
	 set file(phef), file (covariate) from fastlmm_data_ch
	 file(plinks) from  gem_ch_fast
       publishDir "${params.output_dir}/fastlmm", overwrite:true, mode:'copy'
       each this_pheno from ind_pheno_cols_ch
       output:
         file(out)
	 set val(base), val(our_pheno2), file("$out") into fastlmm_manhatten_ch
       script:
	 base = plinks[0].baseName
	 our_pheno = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/-$/,"")
	 our_pheno2 = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/-$/,"").replaceAll(/^[0-9]+-/,"")
	 covar_opt_fast =  (params.covariates) ?  " -covar $covariate" : ""
	 out = "$base-$our_pheno"+".stat"
	 """
	 this_pheno_col=`echo ${this_pheno} | sed 's/-.*//' `
	 $fastlmmc -REML -simType RRM -verboseOut -bfile $base -pheno ${phef} -simLearnType Full -out $out -maxThreads $params.fastlmm_num_cores \
	           $covar_opt_fast  -mpheno \${this_pheno_col} -bfileSim $base
	 """
     }
  }

  // this part is plotting done for any fastlmm mode

  process showFastLmmManhatten {
    publishDir params.output_dir
    input:
      set val(base), val(this_pheno), file(assoc) from fastlmm_manhatten_ch
    output:
      file("${out}*")  into report_fastlmm_ch
    script:
      our_pheno = this_pheno.replaceAll("_","-")
      out = "C051-fastlmm-"+this_pheno
      """
      general_man.py  $assoc $this_pheno ${out} Chromosome Position SNP Pvalue SNPWeight FastLmm
      """
  }


  report_ch = report_ch.flatten().mix(report_fastlmm_ch.flatten())

  // End of FASTLMM
}


/*JT : Case boltlmm => if yes*/



   /*JT Fonction to transforme argument for cofactor in gemma
   @Input 
   args: cofactor args separate by a comma
   infoargs: type of cofactor separate by a comma : 0 for qualitative, 1 for quantitative
   output : cofactor for boltlmm was formating to take account qualitative and quantitative
   */
   def boltlmmCofact(args,infoargs) {
      //Method code
      splargs=args.split(",")
      splinfoargs=infoargs.split(",")
      if(splargs.size() != splinfoargs.size()){
	 System.err.println("args and args type for Boltlmm was not same size : "+args+" "+infoargs)
	 System.exit(-11)
      }
      CofactStr=""
      for (i = 0; i <splargs.size(); i++) {
	  /*0 : for quantitatif */
	  /* 1 for qualitatif*/
	  if     (splinfoargs[i]=='1')  CofactStr +=" --qCovarCol="+splargs[i]
	  else if(splinfoargs[i]=='0')  CofactStr +=" --covarCol="+splargs[i]
	  else{
	     System.err.println("type args for "+splargs[i]+" doesn't know "+ splinfoargs[i]+"\n 1 for quantitatif arguments\n 0 for qualitatif arguments")
	     System.exit(-10)
	  }
      }
      return(CofactStr)
   }

if (params.boltlmm == 1) {

  plink_ch_bolt = Channel.create()
  fam_ch_bolt = Channel.create()
  bim_ch_bolt = Channel.create()
  boltlmm_assoc_ch.separate (plink_ch_bolt, fam_ch_bolt, bim_ch_bolt) { a -> [ a, a[2], a[1]] }
  data_ch_bolt = Channel.fromPath(params.data)
  if (params.covariates)
     covariate_option = "--cov_list ${params.covariates}"
  else
     covariate_option = ""
  process  getBoltPhenosCovar {
    input:
      file(covariates) from data_ch_bolt
      file(fam) from fam_ch_bolt
    output:
      file(phef) into newdata_ch_bolt
      stdout into pheno_cols_ch_bolt
    script:
      base = fam.baseName
      phef = "${base}_fastlmm_n.phe"
      """
      all_covariate.py --data  $covariates --inp_fam  $fam $covariate_option \
                          --pheno ${params.pheno} --phe_out ${phef} --form_out 2
      """
  }

  ind_pheno_cols_ch_bolt = Channel.create()
  check_bolt = Channel.create()
  pheno_cols_ch_bolt.flatMap { list_str -> list_str.split() }.tap ( check_bolt) .set { ind_pheno_cols_ch_bolt }


   if (params.covariates) 
      cov_bolt = boltlmmCofact(params.covariates,params.bolt_covariates_type)
   else
      cov_bolt= ""

   missing_cov=""
   if(params.bolt_use_missing_cov==1)
     missing_cov=" --covarUseMissingIndic "

   pval_head = "P_BOLT_LMM"

   if(params.bolt_ld_score_file){
      ld_score_cmd="--LDscoresFile=$params.bolt_ld_score_file"
      if(params.bolt_ld_scores_col) ld_score_cmd = ld_score_cmd + " --LDscoresCol="+params.bolt_ld_scores_col
   } else 
      ld_score_cmd = "--LDscoresUseChip"

  type_lmm="--lmm"
  
  process doCountNbSnp{
    time   params.big_time
    input :
       file(bim) from bim_ch_bolt
    output :
       stdout into nbsnp_ch_bolt
    script :
      """
      wc -l $bim|awk '{print \$1}'
      """
  }
  /*    nb_snp= CountLinesFile(base+".bim") */
  process doBoltmm{
    cpus params.bolt_num_cores
    memory params.gemma_mem_req
    time   params.big_time
    input:
      set file(plinksbed), file(plinksbim), file(plinksfam) from plink_ch_bolt
      val nb_snp from nbsnp_ch_bolt
      file(phef) from newdata_ch_bolt
    publishDir "${params.output_dir}/boltlmm", overwrite:true, mode:'copy'
    each this_pheno from ind_pheno_cols_ch_bolt
    output:
      file(outReml) into bolt_reml_out
      set val(base), val(our_pheno), file("$out") into bolt_manhatten_ch
    script:
      base = plinksbed.baseName
      our_pheno = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/-$/,"").replaceAll(/^[0-9]+-/,"")
      our_pheno2 = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/-$/,"")
      out     = "$base-$our_pheno2"+".stat"
      outReml = "$base-$our_pheno2"+".reml"
      covar_file_bolt =  (params.covariates) ?  " --covarFile ${phef} " : ""
      model_snp = (nb_snp.toInteger() > 1000000) ? "--modelSnps=.sample.snp" : ""
      """
      shuf -n 950000 $plinksbim | awk '{print \$2}' > .sample.snp 
      bolt $type_lmm --bfile=$base  --phenoFile=${phef} --phenoCol=${our_pheno} --numThreads=$params.bolt_num_cores $cov_bolt $covar_file_bolt --statsFile=$out\
           $ld_score_cmd  $missing_cov --lmmForceNonInf  $model_snp
      bolt  --reml  --bfile=$base  --phenoFile=${phef} --phenoCol=${our_pheno} --numThreads=$params.bolt_num_cores $cov_bolt $covar_file_bolt $missing_cov $model_snp |\
             grep -B 1 -E "^[ ]+h2" > $outReml 
      """
  }

  process showBoltmmManhatten {
    publishDir params.output_dir
    input:
      set val(base), val(this_pheno), file(assoc) from bolt_manhatten_ch
    output:
      file("${out}*")  into report_bolt_ch
    script:
      our_pheno = this_pheno.replaceAll("_","-")
      out = "C052-boltlmmm-"+this_pheno
      """
      general_man.py  $assoc $this_pheno ${out} CHR BP SNP $pval_head BETA BoltLMM
      """
  }
   report_ch = report_ch.flatten().mix(report_bolt_ch.flatten())

}/*JT End of boltlmm*/


  def newNamePheno(Pheno){
      SplP=Pheno.split(',')
      for (i = 0; i <SplP.size(); i++) {
         SplP[i]=(i+1)+"-"+SplP[i]
      }
      return(SplP)
  }



if (params.gemma == 1) {

  rel_ch_gemma = Channel.create()
  gem_ch_gemma = Channel.create()
  gemma_assoc_ch.separate (rel_ch_gemma, gem_ch_gemma) { a -> [a, a] }

  process getGemmaRel {
    cpus params.gemma_num_cores
    memory params.gemma_mem_req
    time params.big_time
    input:
       file plinks from rel_ch_gemma
    output:
       file("output/${base}.*XX.txt") into rel_mat_ch
    script:
       base = plinks[0].baseName
       famfile=base+".fam"
       """
       export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
       cat $famfile |awk '{print \$1"\t"\$2"\t"0.2}' > pheno
       gemma -bfile $base  -gk ${params.gemma_relopt} -o $base -p pheno -n 3
       """
  }

   
  if (params.covariates)
     covariate_option = "--cov_list ${params.covariates}"
  else
     covariate_option = ""

  ind_pheno_cols_ch = newNamePheno(params.pheno)

  process doGemma{
    cpus params.gemma_num_cores
    memory params.gemma_mem_req
    time   params.big_time
    input:
      file(covariates) from data_ch
      file(rel) from rel_mat_ch
      file(plinks) from  gem_ch_gemma
    each this_pheno from ind_pheno_cols_ch
    publishDir params.output_dir, overwrite:true, mode:'copy'
    output: 
      file("gemma/${out}.log.txt")
      set val(newbase), val(our_pheno), file("gemma/${out}.assoc.txt") into gemma_manhatten_ch
    script:
       our_pheno          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/-$/,"")
       our_pheno2         = this_pheno.replaceAll(/^[0-9]+-/,"")
       our_pheno3         = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/-$/,"").replaceAll(/^[0-9]+-/,"")
       data_nomissing     = "pheno-"+our_pheno+".pheno" 
       list_ind_nomissing = "lind-"+our_pheno+".lind"
       rel_matrix         = "newrel-"+our_pheno+".rel"
       base               =  plinks[0].baseName
       inp_fam            =  base+".fam"
       newbase            =  base+"-"+our_pheno
       newfam             =  newbase+".fam"
       gemma_covariate    = "${newbase}.gemma_cov"
       phef               = "${newbase}_n.phe"
       covar_opt_gemma    =  (params.covariates) ?  " -c $gemma_covariate " : ""
       out                = "$base-$our_pheno"
       """
       list_ind_nomissing.py --data $covariates --inp_fam $inp_fam $covariate_option --pheno $our_pheno3 --dataout  $data_nomissing --lindout $list_ind_nomissing
       gemma_relselind.py  --rel $rel --inp_fam $inp_fam --relout $rel_matrix --lind $list_ind_nomissing
       plink --keep-allele-order --bfile $base --keep $list_ind_nomissing --make-bed --out $newbase 
       all_covariate.py --data  $data_nomissing --inp_fam  $newbase".fam" $covariate_option --cov_out $gemma_covariate \
                          --pheno $our_pheno2 --phe_out ${phef} --form_out 1

       export OPENBLAS_NUM_THREADS=${params.gemma_num_cores} 
       gemma -bfile $newbase ${covar_opt_gemma}  -k $rel_matrix -lmm 1  -n 1 -p $phef -o $out -maf 0.0000001 
       mv output gemma
       """
  } 

  process showGemmaManhatten { 
    publishDir params.output_dir
    input:
      set val(base), val(this_pheno), file(assoc) from gemma_manhatten_ch
    output:
      file("${out}*")  into report_gemma_ch
    script:
      our_pheno = this_pheno.replaceAll("_","-")
      out = "C049$this_pheno"
      """
      gemma_man.py  $assoc $this_pheno ${out} 
      """
  }

  report_ch = report_ch.flatten().mix(report_gemma_ch.flatten())
    
} 


    


if (params.chi2+params.fisher+params.logistic+params.linear > 0) {

   process computeTest {
      cpus num_assoc_cores
      time params.big_time
      input:
       set file('cleaned.bed'),file('cleaned.bim'),file('cleaned.fam') from assoc_ch    
       file (phenof) from pheno_ch
      each test from requested_tests
      each pheno_name from pheno_label_ch
      publishDir "${params.output_dir}/${test}", overwrite:true, mode:'copy'
      output:
        set val(test), val(pheno_name), file("${outfname}.*") into out_ch
      script:
       base = "cleaned"
       pheno_name = pheno_name.replaceFirst("/.*","")
       perm = (params.mperm == 0 ? "" : "mperm=${params.mperm}")
       adjust = (params.adjust ? "--adjust" : "")
       outfname = "${pheno_name}"
       if (params.data == "") {
           pheno_cmd = ""
           out = base
       } else {
           pheno_cmd = "--pheno $phenof --pheno-name $pheno_name "
           if (params.covariates) covariate = "--covar ${phenof} --covar-name ${params.covariates} "
           out = pheno
       }
       template "test.sh"
   }


  log_out_ch = Channel.create()
 
  log_out_ch.subscribe { println "Completed plink test ${it[0]}" }
 
  process drawPlinkResults { 
    input:
    set val(test), val(pheno_name), file(results) from out_ch.tap(log_out_ch)
    output:
      set file("${base}*man*png"), file ("${base}*qq*png"), file("C050*tex") into report_plink
    publishDir params.output_dir
    script:
      base="cleaned"
      """
      plinkDraw.py  C050 $base $test ${pheno_name} $gotcovar png
      """
  }

  report_ch = report_ch.mix(report_plink.flatten())
  
}





def getres(x) {
  def  command1 = "$x"
  def  command2 = "head -n 1"
  def proc1 = command1.execute()
  def proc2 = command2.execute()
  def proc = proc1 | proc2
  proc.waitFor()              
  res ="${proc.in.text}"
  return res.trim()
}

nextflowversion =getres("nextflow -v")
if (workflow.repository)
  wflowversion="${workflow.repository} --- ${workflow.revision} [${workflow.commitId}]"
else
  wflowversion="A local copy of the workflow was used"

report_ch = report_ch.mix(report_pca_ch)

process doReport {
  input:
    file(reports) from report_ch.toList()
  publishDir params.output_dir, overwrite:true, mode:'copy'
  output:
    file("${out}.pdf")
  script:
    out = params.output+"-report"
    these_phenos     = params.pheno
    these_covariates = params.covariates
    config = getConfig()
    images = workflow.container
    texf   = "${out}.tex"
    template "make_assoc_report.py"
}





