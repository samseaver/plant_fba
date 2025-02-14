import logging
#import urllib3
import os
import copy
import re
import json
import glob

#For Bokeh figures
#from bokeh.io import output_file, save
#from bokeh.layouts import grid

#For stats
#import pandas as pd
#import scipy as sp
#import scipy.linalg
#import scipy.stats
#import numpy as np

#from installed_clients.DataFileUtilClient import DataFileUtil
#from installed_clients.KBaseReportClient import KBaseReport

#from plant_fba.core.fetch_plantseed_impl import FetchPlantSEEDImpl
#from plant_fba.core.generate_figure_impl import GenerateFigureImpl

class ReconstructAppImpl:
	@staticmethod
	def _validate_params(params, required, optional=set()):
		"""Validates that required parameters are present. Warns if unexpected parameters appear"""
		required = set(required)
		optional = set(optional)
		pkeys = set(params)
		if required - pkeys:
			raise ValueError("Required keys {} not in supplied parameters"
							 .format(", ".join(required - pkeys)))
		defined_param = required | optional
		for param in params:
			if param not in defined_param:
				logging.warning("Unexpected parameter {} supplied".format(param))

	@staticmethod
	def _convert_search_role(role):

		searchrole = role

		#Remove spaces
		searchrole = searchrole.strip()
		searchrole = searchrole.replace(' ','')

		#Make all lowercase
		searchrole = searchrole.lower()

		#Remove EC and parentheses
		searchrole = re.sub(r'\(ec[\d-]+\.[\d-]\.[\d-]\.[\d-]\)', '', searchrole)

		return searchrole

	def _set_objects(self,params):
		self.genome_obj=params['genome']
		self.template_obj=params['template']

	def reconstruct_metabolism(self,input_params):

		#self._validate_params(self.input_params, 
		#					  {'genome_obj',
		#					   'template_obj'})

		# Retrieve Template, and compile indexes of roles and complexes
		# if('template_ws' not in input_params or input_params['template_ws'] == ''):
		#    input_params['template_ws'] = 'NewKBaseModelTemplates'

		# if('template' not in input_params or input_params['template'] == ''):
		#    input_params['template'] = 'PlantModelTemplate'

		# template_ref = input_params['template_ws']+'/'+input_params['template']
		# template_obj = self.dfu.get_objects({'object_refs': [template_ref]})['data'][0]['data']
						
		searchroles_dict = dict()
		roles_dict = dict()
		for role in self.template_obj['roles']:
			searchrole = self._convert_search_role(role['name'])
			searchroles_dict[searchrole]=role['id']
			roles_dict[role['id']]=role

		complex_dict = dict()
		for cpx in self.template_obj['complexes']:
			complex_dict[cpx['id']]=cpx

		abbrev_cpt_dict=dict()
		for cpt in input_params['cpts']:
			abbrev_cpt_dict[input_params['cpts'][cpt]['abbrev']]=cpt

		#Retrieve Genome annotation as dict
		role_cpt_ftr_dict=dict()
		# genome_ref = input_params['input_ws']+'/'+input_params['input_genome']
		# genome_obj = self.dfu.get_objects({'object_refs': [genome_ref]})['data'][0]['data']
		for feature in self.genome_obj['features']:
			if('functions' in feature and len(feature['functions'])>0):
				for function_comment in feature['functions']:

					# print(function_comment)
					#Split for comments and retrieve compartments
					function_cpt_list = function_comment.split("#")
					for i in range(len(function_cpt_list)):
						function_cpt_list[i]=function_cpt_list[i].strip()

					function = function_cpt_list.pop(0)
					# This regex is very old. As I have control over the annotation
					# of plant genomes, I only use the forward slash to separate
					# roles
					# roles = re.split("\s*;\s+|\s+[\@\/]\s+", function)
					roles = re.split("\s+/\s+",function)
					for role in roles:
						
						searchrole = self._convert_search_role(role)

						if(searchrole not in searchroles_dict):
							if(searchrole != 'unannotated'):
								# These are roles that are in PlantSEED but not yet 
								# included in the template
								pass
							continue
							
						role_id = searchroles_dict[searchrole]

						if(role_id not in role_cpt_ftr_dict):
							role_cpt_ftr_dict[role_id]=dict()

						# Defaults to cytosol
						if(len(function_cpt_list)==0):
							function_cpt_list.append('cytosol')

						for cpt in function_cpt_list:
							abbrev_cpt=cpt
							if(cpt not in abbrev_cpt_dict):
								print("No compartmental abbreviation found for "+cpt)
							else:
								abbrev_cpt = abbrev_cpt_dict[cpt]

							if(abbrev_cpt not in role_cpt_ftr_dict[role_id]):
								role_cpt_ftr_dict[role_id][abbrev_cpt]=dict()

							role_cpt_ftr_dict[role_id][abbrev_cpt][feature['id']]=1

		# Default dictionaries for objects needed for a model reaction
		default_mdlcpt_dict = { 'id': 'u0', 'label': 'unknown',
								'pH': 7, 'potential': 0, 'compartmentIndex': 0,
								'compartment_ref': '~//' }

		default_mdlcpd_dict = { 'id': '', 'charge': 0, 'formula': '',
								'name': '', 'compound_ref': '',
								'modelcompartment_ref': '~/modelcompartments/id/u0' }

		default_mdlrxn_dict = { 'id': '', 'direction': '', 'protons': 0,
								'name': '', 'reaction_ref': '', 'probability': 0,
								'modelcompartment_ref': '',
								'modelReactionReagents': [], 'modelReactionProteins': [] }

		# Lookup dictionaries for compartments and compounds, to avoid duplicating them
		mdlcpts_dict = dict()
		mdlcpds_dict = dict()
		
		# Reaction complexes for the generated table
		rxncplxs_dict = dict()

		# 'template_ref' : template_ref, 'genome_ref' : genome_ref, 
		# Create New, but Empty Plant Reconstruction
		new_model_obj = { 'id' : input_params['id'], 'name' : input_params['name'],
				   			'template_ref':'','genome_ref':'',
				   			'type' : "GenomeScale", 'source' : "KBase", 'source_id' : "PlantSEED_v2", 
							'modelreactions' : [], 'modelcompounds' : [], 'modelcompartments' : [], 'biomasses' : [],
							'gapgens' : [], 'gapfillings' : [] }
		
		for template_rxn in self.template_obj['reactions']:
			if(template_rxn['type'] == 'gapfilling'):
				continue

			template_rxn_cpt = template_rxn['templatecompartment_ref'].split('/')[-1]

			proteins_list = list()
			prots_str_list = list()
			# complex_ref and source are optional fields
			default_protein_dict = {'note':template_rxn['type'], 'complex_ref':'', 'modelReactionProteinSubunits':[]}
			for cpx_ref in template_rxn['templatecomplex_refs']:
				cpx_id = cpx_ref.split('/')[-1]
				model_complex_ref = "~/template/complexes/id/"+cpx_id

				new_protein_dict = copy.deepcopy(default_protein_dict)
				new_protein_dict['complex_ref']=model_complex_ref

				complex_present=False
				subunits_list = list()
				default_subunit_dict = {'role':'','triggering':0,'optionalSubunit':0,'note':'','feature_refs':[]}
				matched_role_dict = dict()

				for cpxrole in complex_dict[cpx_id]['complexroles']:
					role_id = cpxrole['templaterole_ref'].split('/')[-1]

					if(role_id in role_cpt_ftr_dict):

						for role_cpt in role_cpt_ftr_dict[role_id]:
							role_cpt_present=False
							if(template_rxn_cpt == role_cpt and cpxrole['triggering'] == 1):
								complex_present=True
								role_cpt_present=True

							if(role_cpt_present == True):
								new_subunit_dict = copy.deepcopy(default_subunit_dict)
								new_subunit_dict['triggering'] = cpxrole['triggering']
								new_subunit_dict['optionalSubunit'] = cpxrole['optional_role']
								new_subunit_dict['role'] = roles_dict[role_id]['name']

								if(len(roles_dict[role_id]['features'])>0):
									new_subunit_dict['note'] = 'Features characterized and annotated'
								else:
									#This never happens as of Fall 2019
									print("Warning: "+roles_dict[role_id]['name']+" is apparently uncharacterized!")
									new_subunit_dict['note'] = 'Features uncharacterized but annotated'
									pass

								for ftr in role_cpt_ftr_dict[role_id][role_cpt]:
									feature_ref = "~/genome/features/id/"+ftr
									new_subunit_dict['feature_refs'].append(feature_ref)
								
								matched_role_dict[role_id]=1
								subunits_list.append(new_subunit_dict)

					if(role_id not in role_cpt_ftr_dict and template_rxn['type'] == 'universal'):
						# This should still be added, with zero features to indicate the universality of the role in plant primary metabolism
						new_subunit_dict = copy.deepcopy(default_subunit_dict)
						new_subunit_dict['triggering'] = cpxrole['triggering']
						new_subunit_dict['optionalSubunit'] = cpxrole['optional_role']
						new_subunit_dict['role'] = roles_dict[role_id]['name']

						# Un-necessary, but explicitly stated
						new_subunit_dict['feature_refs']=[]

						if(len(roles_dict[role_id]['features'])==0):
							new_subunit_dict['note'] = 'Features uncharacterized and unannotated'
						else:
							new_subunit_dict['note'] = "Features characterized but unannotated"
							
							# This includes UniProt annotation and annotation from other genomes
							if(len(roles_dict[role_id]['features']) > 0 and 'Athaliana' in roles_dict[role_id]['features'][0]):
								print("Missing annotation: ",cpx_id,role_id,roles_dict[role_id]['name'],roles_dict[role_id]['features'])

						matched_role_dict[role_id]=1
						subunits_list.append(new_subunit_dict)
						
				if(complex_present == True):
					# Check to see if members of a detected protein complex are missing
					# and add them if so, to round off the complex
					# This will only happen to a complex that is conditional (see above)
					for cpxrole in complex_dict[cpx_id]['complexroles']:
						role_id = cpxrole['templaterole_ref'].split('/')[-1]
						
						if(role_id not in matched_role_dict):
							print("Gapfilling complex: ",cpx_id,roles_dict[role_id])
							new_subunit_dict = copy.deepcopy(default_subunit_dict)
							new_subunit_dict['triggering'] = cpxrole['triggering']
							new_subunit_dict['optionalSubunit'] = cpxrole['optional_role']
							new_subunit_dict['note'] = "Complex-based-gapfilling"
							subunits_list.append(new_subunit_dict)

				if(len(subunits_list)>0):
					new_protein_dict['modelReactionProteinSubunits']=subunits_list

					# Store features and subunits as complex string for table
					subs_str_list=list()
					for subunit in subunits_list:
						ftrs_str_list=list()
						for ftr_ref in subunit['feature_refs']:
							ftr = ftr_ref.split('/')[-1]
							ftrs_str_list.append(ftr)
						ftr_str = "("+", ".join(ftrs_str_list)+")"
						subs_str_list.append(ftr_str)
					sub_str = "["+", ".join(subs_str_list)+"]"
					prots_str_list.append(sub_str)

				proteins_list.append(new_protein_dict)

			prot_str = ", ".join(prots_str_list)

			# This is important, we need to use role-based annotation to determine whether
			# a reaction should even be added to the model
			if(template_rxn['type'] == 'conditional' and len(proteins_list)==0):
				continue

			# If the check passes, then, here, we instantiate the actual reaction that goes into the model
			new_mdlrxn_id = template_rxn['id']+'0'
			new_mdlcpt_id = template_rxn_cpt+'0'
			base_rxn_id = template_rxn['id'].split('_')[0]

			# For table
			rxncplxs_dict[new_mdlrxn_id]=prot_str

			new_mdlrxn_dict = copy.deepcopy(default_mdlrxn_dict)
			new_mdlrxn_dict['id'] = new_mdlrxn_id

			# new_mdlrxn_dict['name'] = MSD_reactions_dict[base_rxn_id]['abbreviation']
			# if(MSD_reactions_dict[base_rxn_id]['abbreviation'] == ""):
			new_mdlrxn_dict['name']=base_rxn_id

			new_mdlrxn_dict['direction'] = template_rxn['direction']
			new_mdlrxn_dict['reaction_ref']='~/template/reactions/id/'+template_rxn['id']
			new_mdlrxn_dict['modelcompartment_ref']='~/modelcompartments/id/'+new_mdlcpt_id

			#Here we check and instantiate a new modelcompartment
			if(new_mdlcpt_id not in mdlcpts_dict):
				new_mdlcpt_dict = copy.deepcopy(default_mdlcpt_dict)
				new_mdlcpt_dict['id']=new_mdlcpt_id
				new_mdlcpt_dict['label']=input_params['cpts'][template_rxn_cpt]['name']
				new_mdlcpt_dict['compartment_ref']='~/template/compartments/id/'+template_rxn_cpt
				mdlcpts_dict[new_mdlcpt_id]=new_mdlcpt_dict

			#Add Proteins as previously determined
			new_mdlrxn_dict['modelReactionProteins']=proteins_list

			#Add Reagents
			for template_rgt in template_rxn['templateReactionReagents']:
				template_rgt_cpd_cpt_id = template_rgt['templatecompcompound_ref'].split('/')[-1]
				(template_rgt_cpd,template_rgt_cpt)=template_rgt_cpd_cpt_id.split('_')

				#Check and add new model compartment 
				new_mdlcpt_id = template_rgt_cpt+'0'
				if(new_mdlcpt_id not in mdlcpts_dict):
					new_mdlcpt_dict = copy.deepcopy(default_mdlcpt_dict)
					new_mdlcpt_dict['id']=new_mdlcpt_id
					new_mdlcpt_dict['label']=input_params['cpts'][template_rxn_cpt]['name']
					new_mdlcpt_dict['compartment_ref']='~/template/compartments/id/'+template_rgt_cpt
					mdlcpts_dict[new_mdlcpt_id]=new_mdlcpt_dict
			   
				#Add new model compounds
				new_mdlcpd_id = template_rgt_cpd_cpt_id+'0'
				base_cpd_id = template_rgt_cpd_cpt_id.split('_')[0]

				if(new_mdlcpd_id not in mdlcpds_dict):
					new_mdlcpd_dict = copy.deepcopy(default_mdlcpd_dict)
					new_mdlcpd_dict['id']=new_mdlcpd_id
					new_mdlcpd_dict['compound_ref']='~/template/compounds/id/'+template_rgt_cpd
					new_mdlcpd_dict['modelcompartment_ref']='~/modelcompartments/id/'+new_mdlcpt_id
					mdlcpds_dict[new_mdlcpd_id]=new_mdlcpd_dict

				new_rgt_dict = {'coefficient' : template_rgt['coefficient'],
								'modelcompound_ref' : '~/modelcompounds/id/'+new_mdlcpd_id}

				new_mdlrxn_dict['modelReactionReagents'].append(new_rgt_dict)

			new_model_obj['modelreactions'].append(new_mdlrxn_dict)

		#Having populated with list of reactions and biomass (to come), then add all compartments and compounds
		for cpt_id in mdlcpts_dict:
			new_model_obj['modelcompartments'].append(mdlcpts_dict[cpt_id])

		#Last, but key modelcompound is the biomass, need to add it explicitly
		biocpd_id = "cpd11416"
		mdlbiocpd_dict = copy.deepcopy(default_mdlcpd_dict)
		mdlbiocpd_dict['id'] = biocpd_id+'_c0'
		mdlbiocpd_dict['name'] = 'Biomass'
		mdlbiocpd_dict['compound_ref'] = "~/template/compounds/id/"+biocpd_id
		mdlbiocpd_dict['modelcompartment_ref'] = "~/modelcompartments/id/c0"
		mdlcpds_dict[mdlbiocpd_dict['id']] = mdlbiocpd_dict

		for cpd_id in mdlcpds_dict:
			new_model_obj['modelcompounds'].append(mdlcpds_dict[cpd_id])

		default_biomass_dict = { 'id': 'bio1', 'name': 'Plant leaf biomass', 'other': 1,
								 'dna': 0, 'rna': 0, 'protein': 0, 'cellwall': 0,
								 'lipid': 0, 'cofactor': 0, 'energy': 0, 'biomasscompounds': [] }

		default_biocpd_dict = { 'modelcompound_ref' : '', 'coefficient' : 0 }

		for template_biomass in self.template_obj['biomasses']:
			new_template_biomass = copy.deepcopy(default_biomass_dict)
			new_template_biomass['id'] = template_biomass['id']
			new_template_biomass['name'] = template_biomass['name']

			for entry in ['dna','rna','protein','cellwall','lipid','cofactor','energy','specialized','other']:
				new_template_biomass[entry] = template_biomass[entry]

			for template_cpd in template_biomass['templateBiomassComponents']:
				new_biocpd_dict = copy.deepcopy(default_biocpd_dict)
				mdlcpd_id = template_cpd['templatecompcompound_ref'].split('/')[-1]+'0'
				if(mdlcpd_id not in mdlcpds_dict):
					print("Template biomass cpd not found in model:",template_cpd)
					continue
				new_biocpd_dict['modelcompound_ref'] = '~/modelcompounds/id/'+mdlcpd_id
				new_biocpd_dict['coefficient'] = template_cpd['coefficient']
				new_template_biomass['biomasscompounds'].append(new_biocpd_dict)
		
			new_model_obj['biomasses'].append(new_template_biomass)

		# print("Saving metabolic reconstruction")
		# model_ws_object = {'type' : 'KBaseFBA.FBAModel', 'name' : input_params['output_fbamodel'],
		# 				   'data' : new_model_obj }

		# if('output_ws' not in input_params or input_params['output_ws'] == ''):
		# 	input_params['output_ws']=input_params['input_ws']

		# ws_id = self.dfu.ws_name_to_id(input_params['output_ws'])
		# saved_model_list=self.dfu.save_objects({'id':ws_id,'objects':[model_ws_object]})[0]
		return new_model_obj

def main():

	# Load test data
	test_data_root = os.path.join("..","..","..","test","data")

	genome_path = os.path.join(test_data_root,"Phytozome_Genomes_Athaliana_TAIR10.Annotated.json")
	genome_fh = open(genome_path)
	genome_obj = json.load(genome_fh)

	template_path = os.path.join(test_data_root,"PlantSEED_Biomass_Template.json")
	template_fh = open(template_path)
	template_obj = json.load(template_fh)

	# Compile compartment information
	compartments = dict()
	with open(os.path.join(test_data_root,'PlantSEED_Compartments.json')) as fh:
		compartments = json.load(fh)

	reconstruct_app = ReconstructAppImpl()
	reconstruct_app._set_objects({'genome':genome_obj,'template':template_obj})

	obj_name='test'
	input_params={'id':obj_name,'name':obj_name,'cpts':compartments}
	metabolism_obj = reconstruct_app.reconstruct_metabolism(input_params)

	msd_biochem_path = "/Users/seaver/Projects/ModelSEEDDatabase/Biochemistry/"
	search_path = os.path.join(msd_biochem_path,"compound_*.json")
	cpds_dict = dict()
	for compounds_file in sorted(glob.glob(search_path)):
		with open(compounds_file) as json_file_handle:
			cpds_list = json.load(json_file_handle)
			for cpd_obj in cpds_list:
				cpds_dict[cpd_obj['id']]=cpd_obj

	for mdlcpd in metabolism_obj['modelcompounds']:
		base_cpd_id = mdlcpd['id'].split('_')[0]
		if(base_cpd_id not in cpds_dict):
			continue

		mdlcpd['name'] = cpds_dict[base_cpd_id]['name']
		mdlcpd['charge'] = float(cpds_dict[base_cpd_id]['charge'])
		mdlcpd['formula'] = cpds_dict[base_cpd_id]['formula']

		if(mdlcpd['formula'] is None):
			mdlcpd['formula'] = ""

		if(' ' in mdlcpd['formula']):
			print(mdlcpd['id'],mdlcpd['formula'])
			
	output_path=os.path.join(test_data_root,obj_name+'.json')
	with open(output_path,'w') as mofh:
		mofh.write(json.dumps(metabolism_obj, indent=4, sort_keys=True))

if(__name__ == "__main__"):
	main()

"""

		# Fetch and parse biochemistry data
		with open(os.path.join("/kb/module/ModelSEEDDatabase",
							   "Biochemistry",
							   "reactions.json")) as msd_rxn_fh:
			MSD_reactions = json.load(msd_rxn_fh)
		MSD_reactions_dict = dict()
		for entry in MSD_reactions:
			MSD_reactions_dict[entry['id']]=entry

		with open(os.path.join("/kb/module/ModelSEEDDatabase",
							   "Biochemistry",
							   "compounds.json")) as msd_rxn_fh:
			MSD_compounds = json.load(msd_rxn_fh)
		MSD_compounds_dict = dict()
		for entry in MSD_compounds:
			MSD_compounds_dict[entry['id']]=entry

		#Compose report string
		html_string="<html><head><title>Reconstruct Plant Metabolism Report</title></head><body>"
		html_string+="<h2>Reconstruct Plant Metabolism Report</h2>"
		html_string+="<p>The \"Reconstruct Plant Metabolism\" app has finished running, "
		html_string+="reconstructing the primary metabolism from the "
		html_string+="enzymatic annotations in "+input_params['input_genome']+"</p>"
		html_string+="<p>Below we present the table of compartmentalized reactions in the metabolic reconstruction, "
		html_string+="it is similar to what you can see in the FBAModel viewer widget that appears "
		html_string+="below the report, but it has some additional information. Each row in the table is unique "
		html_string+="to each combination of reaction and compartment.</p>"
		html_string+="<p><ul>"
		html_string+="<li><b>Subsystems and Classes:</b> The table contains the metabolic subsystems and "
		html_string+="the general class of metabolism they fall into.</li>"
		html_string+="<li><b>Metabolic functions and EC numbers:</b> The table contains the original enzymatic "
		html_string+="annotation ('Roles') and their EC numbers that were associated with each biochemical reaction.</li>"
		html_string+="<li><b>Complexes:</b> The table contains the genes that were annotated with the metabolic functions. "
		html_string+="These genes that are associated with each reaction can be seen in the FBAModel viewer widget, but here "
		html_string+=" one can see how they may be organized into protein complexes. Each set of parentheses '()' "
		html_string+="represents a single protein subunit (which may be the entire enzyme, or part of a large enzymatic "
		html_string+="complex). Each set of square brackets '[]' represents an entire enzyme, regardless of how many "
		html_string+="subunits it consists of. Each reaction may be catalyzed by different enzymes, each in turn composed "
		html_string+="of different subunits. The complexes reflect how the enzymes were curated in <i>Arabidopsis thaliana</i> "
		html_string+=" so if any complex is shown to be empty, this means that the enzymatic annotation was not propagated "
		html_string+="from the original Arabidopsis gene. The original Arabidopsis curation also included protein localization "
		html_string+="so if a reaction has empty complexes in some compartments as opposed to others, this is an indication "
		html_string+="that annotation was only propagated for some localized Arabidopsis enzymes, and not others."
		html_string+="</ul></p>"

		# Fetch PlantSEED Data
		with open(os.path.join("/kb/module/PlantSEED",
							   "Data/PlantSEED_v3",
							   "PlantSEED_Roles.json")) as plsd_fh:
			PS_Roles = json.load(plsd_fh)

		plantseed = FetchPlantSEEDImpl()
		reactions_data = plantseed.fetch_reactions(PS_Roles)

		table = GenerateTableImpl()
		table_html_string = table.generate_table(reactions_data, complexes=rxncplxs_dict)

		with open(os.path.join('/kb/module/data','app_report_templates','integrate_abundances_report_tables_template.html')) as report_template_file:
			report_template_string = report_template_file.read()

		# Generate and insert html Title
		report_template_string = report_template_string.replace('*TITLE*', input_params['output_fbamodel'])

		# Insert html table
		table_report_string = report_template_string.replace('*TABLES*', html_string+table_html_string)

		#Make folder for report files
		uuid_string = str(uuid.uuid4())
		report_file_path=os.path.join(self.shared_folder,uuid_string)
		os.mkdir(report_file_path)

		#Write html files
		with open(os.path.join(report_file_path,"index.html"),'w') as index_file:
			index_file.write(table_report_string)

		#Cache it in shock as an archive
		upload_info = self.dfu.file_to_shock({'file_path': report_file_path,
											  'pack': 'zip'})

		#Prepare report parameters
		report_params = { 'direct_html_link_index' : 0, #Use to refer to index of 'html_links'
						  'workspace_name' : input_params['input_ws'],
						  'report_object_name' : 'plant_fba_' + uuid_string,
						  'objects_created' : [],
						  'html_links' : [] }

		#Html Link object
		html_link = {'shock_id' : upload_info['shock_id'],
					 'name' : 'index.html',
					 'label' : 'html files',
					 'description' : 'HTML files'}
		report_params['html_links'].append(html_link)

		#Objects created object
		saved_model_ref = "{}/{}/{}".format(saved_model_list[6],saved_model_list[0],saved_model_list[4])
		saved_model_desc = "FBAModel: "+input_params['output_fbamodel']
		report_params['objects_created'].append({'ref':saved_model_ref,'description':saved_model_desc})

		kbase_report_client = KBaseReport(self.callback_url, token=self.token)
		report_client_output = kbase_report_client.create_extended_report(report_params)

		output_report=dict()
		output_report['report_name']=report_client_output['name']
		output_report['report_ref']=report_client_output['ref']            
"""
