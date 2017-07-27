#!/usr/bin/env bash

tax_id=$1
output_file=$2
format=$3

# Names and taxononmy
entry_name="entry name"
gene_names="genes"
gene_names_primary="genes(PREFERRED)"
gene_names_synonym="genes(ALTERNATIVE)"
gene_names_ordered_locus="genes(OLN)"
gene_names_ORF="genes(ORF)"
organism="organism"
organism_ID="organism-id"
protein_names="protein names"
proteomes="proteome"
taxonomic_lineage="lineage(ALL)"
virus_hosts="virus hosts"

# Sequences
fragment="fragment"
gene_encoded_by="encodedon"
alternative_products="comment(ALTERNATIVE PRODUCTS)"
erroneous_gene_model_prediction="comment(ERRONEOUS GENE MODEL PREDICTION)"
erroneous_initiation="comment(ERRONEOUS INITIATION)"
erroneous_termination="comment(ERRONEOUS TERMINATION)"
erroneous_translation="comment(ERRONEOUS TRANSLATION)"
frameshift="comment(FRAMESHIFT)"
mass_spectrometry="comment(MASS SPECTROMETRY)"
polymorphism="comment(POLYMORPHISM)"
RNA_editing="comment(RNA EDITING)"
sequence_caution="comment(SEQUENCE CAUTION)"
length="length"
mass="mass"
sequence="sequence"
alternative_sequence="feature(ALTERNATIVE SEQUENCE)"
natural_variant="feature(NATURAL VARIANT)"
non_adjacent_residues="feature(NON ADJACENT RESIDUES)"
Non_standard_residue="feature(NON STANDARD RESIDUE)"
non_terminal_residue="feature(NON TERMINAL RESIDUE)"
sequence_conflict="feature(SEQUENCE CONFLICT)"
sequence_uncertainty="feature(SEQUENCE UNCERTAINTY)"
sequence_version="version(sequence)"


# Function
EC_number="ec"
absorption="comment(ABSORPTION)"
catalytic_activity="comment(CATALYTIC ACTIVITY)"
cofactor="comment(COFACTOR)"
enzyme_regulation="comment(ENZYME REGULATION)"
function_CC="comment(FUNCTION)"
kinetics="comment(KINETICS)"
pathway="comment(PATHWAY)"
redox_potential="comment(REDOX POTENTIAL)"
temperature_dependence="comment(TEMPERATURE DEPENDENCE)"
pH_dependence="comment(PH DEPENDENCE)"
active_site="feature(ACTIVE SITE)"
binding_site="feature(BINDING SITE)"
DNA_binding="feature(DNA BINDING)"
metal_binding="feature(METAL BINDING)"
nucleotide_binding="feature(NP BIND)"
site="feature(SITE)"

# Miscellaneous
annotation_score="annotation score"
features="features"
caution="comment(CAUTION)"
miscellaneous_CC="comment(GENERAL)"
keywords="keywords"
matched_text="context"
protein_existence="existence"
tools="tools"
reviewed="reviewed"

# Interaction
subunit_structure_CC="comment(SUBUNIT)"
interacts_with="interactor"

# Expression
developmental_stage="comment(DEVELOPMENTAL STAGE)"
induction="comment(INDUCTION)"
tissue_specificity="comment(TISSUE SPECIFICITY)"

# Gene Ontology
gene_ontology_GO="go"
gene_ontology_biological_process="go(biological process)"
gene_ontology_molecular_function="go(molecular function)"
gene_ontology_cellular_component="go(cellular component)"
gene_ontology_IDs="go-id"

# Pathology and biotech
allergenic_properties="comment(ALLERGEN)"
biotechnological_use="comment(BIOTECHNOLOGY)"
disruption_phenotype="comment(DISRUPTION PHENOTYPE)"
involvement_in_disease="comment(DISEASE)"
pharmaceutical_use="comment(PHARMACEUTICAL)"
toxic_dose="comment(TOXIC DOSE)"

# Subcellular location
subcellular_location_CC="comment(SUBCELLULAR LOCATION)"
intramembrane="feature(INTRAMEMBRANE)"
topological_domain="feature(TOPOLOGICAL DOMAIN)"
transmembrane="feature(TRANSMEMBRANE)"

# PTM / Processing
post_translational_modification="comment(PTM)"
chain="feature(CHAIN)"
cross_link="feature(CROSS LINK)"
disulfide_bond="feature(DISULFIDE BOND)"
glycosylation="feature(GLYCOSYLATION)"
initiator_methionine="feature(INITIATOR METHIONINE)"
lipidation="feature(LIPIDATION)"
modified_residue="feature(MODIFIED RESIDUE)"
peptide="feature(PEPTIDE)"
propeptide="feature(PROPEPTIDE)"
signal_peptide="feature(SIGNAL)"
transit_peptide="feature(TRANSIT)"

# Structure
_3D="3d"
beta_strand="feature(BETA STRAND)"
helix="feature(HELIX)"
turn="feature(TURN)"

# Publications
mapped_PubMed_ID="citationmapping"
PubMed_ID="citation"

# Date of
date_of_creation="created"
date_of_last_modification="last-modified"
date_of_last_sequence_modification="sequence-modified"
entry_version="version(entry)"

# Family and domains
domain_CC="comment(DOMAIN)"
sequence_similarities="comment(SIMILARITY)"
protein_families="families"
coiled_coil="feature(COILED COIL)"
compositional_bias="feature(COMPOSITIONAL BIAS)"
domain_FT="feature(DOMAIN EXTENT)"
motif="feature(MOTIF)"
region="feature(REGION)"
repeat="feature(REPEAT)"
zinc_finger="feature(ZINC FINGER)"

# Taxononomic lineage
taxonomic_lineage_all="lineage(all)"

# Taxonomic identifier
taxonomic_identifier_all="lineage-id(all)"

query_list=("${entry_name}" \
"${gene_names}" \
"${gene_names_primary}" \
"${gene_names_synonym}" \
"${gene_names_ordered_locus}" \
"${gene_names_ORF}" \
"${organism}" \
"${organism_ID}" \
"${protein_names}" \
"${proteomes}" \
"${taxonomic_lineage}" \
"${virus_hosts}" \
"${fragment}" \
"${gene_encoded_by}" \
"${alternative_products}" \
"${erroneous_gene_model_prediction}" \
"${erroneous_initiation}" \
"${erroneous_termination}" \
"${erroneous_translation}" \
"${frameshift}" \
"${mass_spectrometry}" \
"${polymorphism}" \
"${RNA_editing}" \
"${sequence_caution}" \
"${length}" \
"${mass}" \
"${sequence}" \
"${alternative_sequence}" \
"${natural_variant}" \
"${non_adjacent_residues}" \
"${Non_standard_residue}" \
"${non_terminal_residue}" \
"${sequence_conflict}" \
"${sequence_uncertainty}" \
"${sequence_version}" \
"${EC_number}" \
"${absorption}" \
"${catalytic_activity}" \
"${cofactor}" \
"${enzyme_regulation}" \
"${function_CC}" \
"${kinetics}" \
"${pathway}" \
"${redox_potential}" \
"${temperature_dependence}" \
"${pH_dependence}" \
"${active_site}" \
"${binding_site}" \
"${DNA_binding}" \
"${metal_binding}" \
"${nucleotide_binding}" \
"${site}" \
"${annotation_score}" \
"${features}" \
"${caution}" \
"${miscellaneous_CC}" \
"${keywords}" \
"${matched_text}" \
"${protein_existence}" \
"${tools}" \
"${reviewed}" \
"${subunit_structure_CC}" \
"${interacts_with}" \
"${developmental_stage}" \
"${induction}" \
"${tissue_specificity}" \
"${gene_ontology_GO}" \
"${gene_ontology_biological_process}" \
"${gene_ontology_molecular_function}" \
"${gene_ontology_cellular_component}" \
"${gene_ontology_IDs}" \
"${allergenic_properties}" \
"${biotechnological_use}"\
"${disruption_phenotype}" \
"${involvement_in_disease}" \
"${pharmaceutical_use}" \
"${toxic_dose}" \
"${subcellular_location_CC}" \
"${intramembrane}" \
"${topological_domain}" \
"${transmembrane}" \
"${post_translational_modification}" \
"${chain}" \
"${cross_link}" \
"${disulfide_bond}" \
"${glycosylation}" \
"${initiator_methionine}" \
"${lipidation}" \
"${modified_residue}" \
"${peptide}" \
"${propeptide}" \
"${signal_peptide}" \
"${transit_peptide}" \
"${_3D}" \
"${beta_strand}" \
"${helix}" \
"${turn}" \
"${mapped_PubMed_ID}" \
"${PubMed_ID}" \
"${date_of_creation}" \
"${date_of_last_modification}" \
"${date_of_last_sequence_modification}" \
"${entry_version}" \
"${domain_CC}" \
"${sequence_similarities}" \
"${protein_families}" \
"${coiled_coil}" \
"${compositional_bias}" \
"${domain_FT}" \
"${motif}" \
"${region}" \
"${repeat}" \
"${zinc_finger}" )

query_string=$( echo $( IFS=$','; echo "${query_list[*]}" ) )

if [ ${format} -eq fasta ]; then
    query_sting=$sequence;
fi

url="http://www.uniprot.org/uniprot/?query=organism:${tax_id}&sort=score&columns=${query_string}&format=${format}"

wget -O $output_file "$url" 


