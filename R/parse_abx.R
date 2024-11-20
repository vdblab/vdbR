# Nov/12/2021
#
# Unify the way we classify drug exposures

# as of 2021-11-12, this was all the drugs that were present in the most recent pull of the tasks table.
drug_classes <- read.table(header = TRUE, sep = ",", text = "drug_name_clean,class,spectrum,betalactam
ciprofloxacin,quinolone,broad,0
piperacillin_tazobactam,penicillins,broad,1
vancomycin_oral,glycopeptide_oral,grampos,0
vancomycin_iv,glycopeptide_iv,grampos,0
vancomycin_unknown,glycopeptide_unknown,grampos,0
vancomycin,glycopeptide,grampos,0
imipenem_cilastatin,carbapenem,broad,1
cefepime,cephem,broad,1
metronidazole,nitroimidazole,broad,0
meropenem,carbapenem,broad,1
cefazolin,cephem,broad,1
sulfamethoxazole_trimethoprim,sulfonamide,broad,0
daptomycin,lipopeptide,grampos,0
linezolid,oxazolidone,grampos,0
levofloxacin,quinolone,broad,0
azithromycin,macrolide,broad,0
erythromycin,macrolide,broad,0
bacitracin,polypeptide,grampos,0
aztreonam,monobactam,gramneg,1
cefadroxil,cephem,broad,1
gentamicin,aminoglycoside,broad,0
tigecycline,tetracycline,broad,0
polymixin b,polymyxin,gramneg,0
ampicillin,penicillins,broad,1
amikacin,aminoglycoside,moderate,0
amoxicillin,penicillins,moderate,1
ceftriaxone,cephem,broad,1
cefixime,cephem,broad,1
ceftazidime,cephem,broad,1
clindamycin,lincosamides,broad,0
cephalexin,cephem,broad,1
penicillin,penicillins,broad,0
doxycycline,tetracycline,broad,0
ertapenem,carbapenem,broad,1
dapsone,sulfone,other,0
cefuroxime,cephem,broad,1
acyclovir,nucleoside,antiviral,0
fluconazole,azole,antifungal,0
amphotericin_b_liposome,macrolide,antifungal,0
micafungin,polyene,antifungal,0
voriconazole,azole,antifungal,0
levofloxacin_(levaquin),quinolone,broad,0
posaconazole,azole,antifungal,0
ticarcillin_clavulanate,penicillins,broad,1
atovaquone,quinone,antifungal,0
famciclovir,nucleoside,antiviral,0
oseltamivir,neuraminidase,antiviral,0
amoxicillin_clavulanate,penicillins,moderate,1
cefotetan,cephem,broad,1
cidofovir,DNApol_inhib,antiviral,0
valacyclovir,DNApol_inhib,antiviral,0
entecavir,nucleoside,antiviral,0
foscarnet,DNApol_inhib,antiviral,0
valganciclovir,nucleoside,antiviral,0
isoniazid,prodrug,narrow,0
rifampin,ansamycin,broad,0
pyrazinamide,prodrug,narrow,0
ethambutol,cellwallinhibitor,narrow,0
ganciclovir,nucleoside,antiviral,0
ketoconazole,azole,antifungal,0
ampicillin_sulbactam,penicillins,broad,1
ivermectin,calciumchannelblock,antiparasitic,0
tenofovir,nucleoside,antiviral,0
oxacillin,penicillins,narrow,1
polymyxin_b_sulfate,polymyxin,gramneg,0
pyrimethamine,dhfr_inhib,antiparasitic,0
micafungin_invest,polyene,antifungal,0
primaquine,quinoline,antifungal,0
lamivudine,nucleoside,antiviral,0
rifabutin,ansamycin,broad,0
sulfadiazine,folate_synth_inhib,broad,0
nitrofurantoin,nitrofuran,broad,0
tigecycline_invest,tetracycline,broad,0
rifaximin,RNApol_inhib,broad,0
dicloxacillin,penicillins,narrow,1
demeclocycline,tetracycline,broad,0
ritonavir,protease_inhib,antiviral,0
emtricitabine_tenofovir,nucleoside,antiviral,0
maraviroc_invest,entry_inhib,antiviral,0
darunavir,protease_inhib,antiviral,0
isavuconazonium,azole,antifungal,0
cefpodoxime,cephem,broad,1
paromomycin,protein_synth_inhib,broad,0
isavuconazonium_sulfate_invest,azole,antifungal,0
ceftaroline,cephem,broad,1
nitazoxanide,thiazolides,antiviral,0
cefoxitin,cephem,broad,1
letermovir,terminase_inhib,antiviral,0
minocycline,tetracycline,broad,0
meropenem_vaborbactam,carbapenem,broad,1
letermovir_(mk_8228)_invest,terminase_inhib,antiviral,0
flucytosine,DNA_synthesis_inhib,antifungal,0
rimantadine,DNApol_inhib,antiviral,0
amphotericin_b,macrolide,antifungal,0
neomycin,aminoglycoside,broad,0
azithromycin_ophthalmic,macrolide,broad,0
abacavir_dolutegravir_lamivudine,nucleoside,antiviral,0
abacavir,nucleoside,antiviral,0
dolutegravir,DNAstrand_inhib,antiviral,0
letermovir_(mk_8228)_or_placebo_invest,terminase_inhib,antiviral,0
")
#' get_abx_from_tasks
#' @name get_abx_from_tasks
#' @param tasks_path path to tasks pull in rds format, such as /deep_sequencing/Clinical Annotation/Abx_data_TsoniIDB_pull_2020_07-06/allo_tasks_2020-07-10.rds
#' @export
#'
get_abx_from_tasks <- function(tasks_path) {
  allo_tasks_raw <- readRDS(tasks_path)

  # pull all the antibiotics, antifungals, antiparasitics, and antivirals under the anti-infective category
  # clean up the name by swapping  spaces, dashes for underscores, and make it all lowercase
  antiinf_allo_tasks <- allo_tasks_raw %>%
    dplyr::filter("anti-infectives" == med_class1) %>%
    dplyr::mutate(drug_name_clean = gsub(" |-", "_", tolower(drug_name))) %>%
    dplyr::left_join(., drug_classes, by = "drug_name_clean") %>%
    dplyr::rename(mrn = MRN) %>%
    dplyr::filter(grepl("Performed", task_status))

  empty_class <- antiinf_allo_tasks %>%
    dplyr::select(mrn, drug_name_clean, class, spectrum) %>%
    dplyr::filter(is.na(class))

  if (nrow(empty_class) > 0) {
    warning("some anti-infectives not categorized; please submit an issue/PR with the following info:")
    print(empty_class)
  }
  return(antiinf_allo_tasks)
}

