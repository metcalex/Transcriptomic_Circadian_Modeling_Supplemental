import json, sys, warnings, os, math
import matplotlib
import matplotlib.pyplot as plt
import cobra
import importlib
import copy
from pprint import pprint

from scipy.interpolate import interp1d

import re
import time

import pandas as pd

import numpy as np

import signal

import stopit

#Sets an handler function, you can comment it if you don't need it.
class CustomTimeExit(Exception):
	pass

def handler_function(signum, frame):
	print('Ended')
	raise CustomTimeExit('Ended')

signal.signal(signal.SIGALRM,handler_function)

class timeout:
    def __init__(self, seconds=2, error_message='Timeout'):
        self.seconds = seconds
        self.error_message = error_message
    def handle_timeout(self, signum, frame):
        raise TimeoutError(self.error_message)
    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.alarm(self.seconds)
    def __exit__(self, type, value, traceback):
        signal.alarm(0)

import sys
sys.path.append('./')

import AICcModels
importlib.reload(AICcModels)

def make_storage_directory(storage_directory):
	if not os.path.isdir(storage_directory):
		os.mkdir(storage_directory)

def model_chunker(model_name, run_config_dict):
	#first check to see if model exists, if it does, return it
	if os.path.exists('Outputs/{}'.format(run_config_dict["Chunked Model"])):
		print('Loading previously chunked model at Outputs/{}'.format(run_config_dict["Chunked Model"]))
		with warnings.catch_warnings():
			warnings.simplefilter("ignore")
			return cobra.io.read_sbml_model('Outputs/{}'.format(run_config_dict["Chunked Model"]))

	#chunk model from information
	print('Creating new chunked model from {} ...'.format(model_name), end='\t')
	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
		model=cobra.io.read_sbml_model(model_name)
		cobra.io.write_sbml_model(model,'Outputs/Orig_Model.xml')
	chunk_def_dict=json.load(open('chunk_definitions.json'))
	mass_dict=json.load(open('MetaboliteMassList.json'))
	output_model=copy.deepcopy(model)
	for reaction_id in chunk_def_dict:
		mass=0
		for metabolite in chunk_def_dict[reaction_id]['Metabolites']:
			mass+=chunk_def_dict[reaction_id]['Metabolites'][metabolite]* float(mass_dict.get(metabolite,0))
			#Add dumping of unprocessessable metabolites
			if float(mass_dict.get(metabolite,0))>0:
				new_export_id='{}_chunk_export'.format(metabolite)
				try:
					output_model.reactions.get_by_id(new_export_id)
				except KeyError:
					new_metabolite_waste_reaction=cobra.Reaction(new_export_id)
					output_model.add_reaction(new_metabolite_waste_reaction)
					new_metabolite_waste_reaction.add_metabolites({metabolite:-1})
					new_metabolite_waste_reaction.__imul__(1000/float(mass_dict.get(metabolite,0)))
					new_metabolite_waste_reaction.name='{} chunk dump'.format(metabolite)
					new_metabolite_waste_reaction.bounds=(0,0)
			else:
				print('{} is massless, cannot be dumped'.format(metabolite))


		print((reaction_id,mass))
		print()
		new_chunk_metabolite=cobra.Metabolite('{}_c'.format(reaction_id),  compartment='c')
		output_model.add_metabolites([new_chunk_metabolite])
		new_chunk_metabolite.name=reaction_id
		print(new_chunk_metabolite.id)
		new_chunk_reaction=cobra.Reaction('{}_formation'.format(reaction_id))
		output_model.add_reaction(new_chunk_reaction)
		new_chunk_reaction.add_metabolites(chunk_def_dict[reaction_id]['Metabolites'])
		new_chunk_reaction.name=chunk_def_dict[reaction_id]['Name']
		new_chunk_reaction.__imul__(-1000/mass)
		new_chunk_reaction.add_metabolites({new_chunk_metabolite:1},combine=False)
		if reaction_id == 'Energy':
			new_chunk_reaction.bounds=(0,1000) #note: this will be adjusted
		else:
			new_chunk_reaction.bounds=(-1000,1000)
			new_chunk_export=cobra.Reaction('{}_export'.format(reaction_id))
			output_model.add_reaction(new_chunk_export)
			new_chunk_export.add_metabolites({new_chunk_metabolite:-1})
			new_chunk_export.bounds=(0,1000)
			new_chunk_reaction.name='{} Export'.format(chunk_def_dict[reaction_id]['Name'])
	cobra.io.write_sbml_model(output_model,'Outputs/{}'.format(run_config_dict["Chunked Model"]))
	print('Saving model to Outputs/{}'.format(run_config_dict["Chunked Model"]))
	return output_model

def gene_rule_dictionary_maker(gene_rule):
	GRRdebugFlag=False
	#print(gene_rule)
	#if 'Cre10.g420350.t1.2' in gene_rule and 'Cre05.g238332.t1.1' in gene_rule and True:
	if 'Cre02.g099850.t1.1' in gene_rule and 'Cre09.g386758.t1.1' in gene_rule and False:
		print('\n\nDEBUG RUN\n\n'+('*'*50)+'\n\n')
		print(gene_rule)
		GRRdebugFlag=True
	#print('\n\nRunning {}'.format(gene_rule))
	#remove leading and lagging spaces
	gene_rule=gene_rule.strip()
	#print(gene_rule)
	#ensure all parantheses have spaces around them
	#match ' *\( *'
	def min_match_spacer(matchobj):
		#print(matchobj)
		#sys.exit()
		#return ' m '
		#print(' {} '.format(matchobj.group(0).strip().replace('(','[').replace(')',']')))
		#sys.exit()
		return ' {} '.format(matchobj.group(0).strip().replace('(','[').replace(')',']'))
	if GRRdebugFlag:
		print()
		print(gene_rule)
	while 'min(' in gene_rule:
		gene_rule=re.sub(r' *min\((?!.*min\().*?\) *',min_match_spacer,gene_rule) #matches min(
		#print(gene_rule)
	if GRRdebugFlag:
		print(gene_rule)
		#sys.exit()
	gene_rule=re.sub(r' *(?<!min)\( *',' ( ',gene_rule) #matches ( but not min(
	#print(gene_rule)
	gene_rule=re.sub(r' *\) *',' ) ',gene_rule)
	#gene_rule=re.sub(r' *(?<!min\(.*?)\) *',' ) ',gene_rule) #matches ) but not the ) in min()
	gene_rule=gene_rule.replace('[','(').replace(']',')')


	#print(gene_rule)

	#generate a word list
	word_list=gene_rule.split(' ')
	if 'min' in gene_rule and 'and' not in gene_rule and GRRdebugFlag:
		pass
		print(word_list)
		#sys.exit()
	#pprint(word_list)
	operation_dictionary={}
	#repeatedly check the word list
	if GRRdebugFlag:
		word_eval_count=0
	while True:
		if GRRdebugFlag:
			print(word_eval_count,end='\t')
			print(word_list)
			word_eval_count+=1
		index=0
		r_paran_loc=-1
		and_loc=-1
		or_loc=-1
		l_paran_loc=-1
		#move along the word list
		while index<len(word_list):
			word=word_list[index]
			if word =='(':
				l_paran_loc=index
			#elif word =='min(':
			#	l_paran_loc=index
			elif word.lower() == 'and':
				and_loc=index
			elif word.lower()== 'or':
				or_loc=index
			elif word == ')':
				r_paran_loc=index
			if r_paran_loc>-1 or (l_paran_loc==-1 and (and_loc>-1 or or_loc>-1) and word_list[index+1] != '('):
				#print('Here')
				#print(word_list[index])
				#print((r_paran_loc>-1, (l_paran_loc==-1 and (and_loc>-1 or or_loc>-1) and word_list[index+1] != '(')))
				break
			index+=1
		#build the new word list
		if index==len(word_list):
			#print('Returning {}'.format(''.join(word_list)))
			if GRRdebugFlag:
				print(''.join(word_list))
				#sys.exit()
			return ''.join(word_list)
		new_word_list=[]
		if r_paran_loc==-1:
			if and_loc>-1:
				new_word_list.extend(word_list[:and_loc-1])
				new_word_list.append('min({},{})'.format(word_list[and_loc-1], word_list[and_loc+1]))
				if len(word_list)>and_loc+1:
					new_word_list.extend(word_list[and_loc+2:])
			elif or_loc>-1:
				new_word_list.extend(word_list[:or_loc-1])
				new_word_list.append('({}+{})'.format(word_list[or_loc-1], word_list[or_loc+1]))
				if len(word_list)>or_loc+1:
					new_word_list.extend(word_list[or_loc+2:])
		else:
			#print((word_list[l_paran_loc:r_paran_loc]))
			new_word_list.extend(word_list[:l_paran_loc])
			new_word_list.append(gene_rule_dictionary_maker(' '.join(word_list[l_paran_loc+1:r_paran_loc])))
			new_word_list.extend(word_list[r_paran_loc+1:])
		#print(new_word_list)
		word_list=new_word_list
		#time.sleep(0.5)


def gene_value_finder(model_list, gene_info, current_time):
	if max(gene_info['y'])>1:
		return float(model_dict[gene_info['Model']['BestFunc']](current_time,*gene_info['Model']['popt']))
	else:
		return 1

def reaction_limit_evaluator(reaction_gene_dict, gene_value_dict):
	#this evaluates the reaction limits for a known gene set
	output_dict={}
	for rxn_id in reaction_gene_dict:
		rule_string=reaction_gene_dict[rxn_id]['Rule']
		print(rule_string)
		for gene in reaction_gene_dict[rxn_id]['Genes']:
			rule_string=rule_string.replace(gene,str(gene_value_dict[gene]))
		print(rule_string)
		#print()
		output_dict[rxn_id]=eval(rule_string)
	return output_dict


def gene_value_dict_maker(current_time, total_gene_set, conversion_file, gene_max_dict, model_parameter_file, run_config_dict):
	#This function caculates the value for all genes
	gene_dict={}
	for gene in total_gene_set:
		gene_knockout_flag=False
		#print(gene)
		for knockout in run_config_dict['Knockout Gene List']:
			if knockout in gene:
				gene_knockout_flag=True
				#print('{} is knocking out {}')
				#print(gene)
				#print(knockout)
				#sys.exit()
		if not gene_knockout_flag:
			if gene in conversion_file:
				gene_key_list=conversion_file[gene]
				min_gene=np.inf
				#checks for multiple genes
				for gene_key in gene_key_list:
					gene_value=gene_value_finder(model_dict, model_parameter_file[gene_key],current_time)
					min_gene=min(min_gene,gene_value)
				gene_dict[gene]=min_gene/gene_max_dict[gene]

			elif '.' in gene:
				gene_key='.'.join(gene.split('.')[:-2])
				gene_value=gene_value_finder(model_dict, model_parameter_file[gene_key],current_time)
				gene_dict[gene]=gene_value/gene_max_dict[gene]
		else:
			gene_dict[gene]=0

	return gene_dict

def model_bounds_adjuster(chunked_model, fixed_chunked_model,gene_value_dict, reaction_gene_dict, max_flux, current_time, run_config_dict):
	reaction_exclusion_list=['EX_photonVis_e']
	reaction_exclusion_list.extend(list(chunked_model.medium.keys()))
	for reaction in chunked_model.reactions:
		if reaction.id not in reaction_exclusion_list:
			#print(reaction.id)

			bnds=fixed_chunked_model.reactions.get_by_id(reaction.id).bounds
			if "Use Transcript Boolean" in run_config_dict and run_config_dict["Use Transcript Boolean"]:
				reaction.bounds=bnds
			else:
				rule_string=reaction_gene_dict[reaction.id]['Rule']
				#print(rule_string)
				for gene in reaction_gene_dict[reaction.id]['Genes']:
					rule_string=rule_string.replace(gene,str(gene_value_dict.get(gene,1)))
				if rule_string!='':
					ret_val=min(max_flux*eval(rule_string),max_flux)
				else:
					ret_val=max_flux
				ret_val=max(ret_val,0)
				#new_bounds=[]
				if bnds[0]>0 or bnds[1]<0 or bnds==(0,0):
					new_bounds=bnds
				elif bnds[0]<0:
					new_bounds=(-ret_val,ret_val)
				else:
					new_bounds=(0,ret_val)
				reaction.bounds=new_bounds
	return chunked_model
		#sys.exit()

def time_stringer(input_time):
	return '{0:2.5f}'.format(float(input_time))

def exchange_dict_comparer(ex_dict, time_feas_info):
	media_non_photon_check_flag=True
	for media_key in ex_dict:
		if media_key != 'EX_photonVis_e':
			if media_key not in time_feas_info['Exchange Bounds']:
				media_non_photon_check_flag=False
				print('{} not in time_feas dict'.format(media_key))
			elif ex_dict[media_key]!=time_feas_info['Exchange Bounds'][media_key]:
				check_lb_flag=time_feas_info['Exchange Bounds'][media_key][0]!= ex_dict[media_key][0] and (time_feas_info['Exchange Bounds'][media_key][0]>=1e-5 or ex_dict[media_key][0]>=1e-5)
				check_ub_flag=time_feas_info['Exchange Bounds'][media_key][1]!= ex_dict[media_key][1] and (time_feas_info['Exchange Bounds'][media_key][1]<=1e-5 or ex_dict[media_key][1]<=1e-5)
				if check_lb_flag or check_ub_flag:
					media_non_photon_check_flag=False
					print('Major difference between {} values'.format(media_key))
					print(ex_dict[media_key])
					print(time_feas_info['Exchange Bounds'].get(media_key, None))
	return media_non_photon_check_flag

def model_feasibility_tester(chunked_model,fixed_chunked_model, total_gene_set, conversion_file, gene_max_dict, reaction_gene_dict, chunk_feasibilty_dict, run_config_dict, base_chunk_feasibility_dict, time_set=np.linspace(0,24,481)):
	#This function checks what is required to build the cell's chunks at different times.
	met_list=[]
	for reaction in chunked_model.reactions:
		if 'CH_' in reaction.id and '_export' in reaction.id:
			met_list.append('CH_{}'.format(reaction.id.split('_')[1]))
	#pprint(chunk_feasibilty_dict)
	#chunk_feasibilty_dict={}
	pprint(met_list)

	current_time=time_set[0]
	#set model to first time point
	gene_value_dict=gene_value_dict_maker(float(current_time), total_gene_set, conversion_file, gene_max_dict, model_parameter_file, run_config_dict)
	chunked_model=model_bounds_adjuster(chunked_model,fixed_chunked_model, gene_value_dict, reaction_gene_dict, max_flux, float(current_time), run_config_dict)

	#check each time
	overall_run_flag=False
	for current_time in time_set:

		string_time=time_stringer(current_time)
		chunked_model, light_level=light_and_media_adjuster(run_config_dict, chunked_model, current_time)
		#make the exchange dict
		ex_dict={}
		for ex_rxn in chunked_model.exchanges:
			ex_dict[ex_rxn.id]=list(ex_rxn.bounds)

		#check if all the exchanges match
		if string_time in chunk_feasibilty_dict:
			time_feas_info=chunk_feasibilty_dict[string_time]
			media_non_photon_check_flag=exchange_dict_comparer(ex_dict, time_feas_info)
		else:
			media_non_photon_check_flag=False
		print(media_non_photon_check_flag,end='\t')

		run_flag=True
		if string_time in chunk_feasibilty_dict and media_non_photon_check_flag:
			run_info=chunk_feasibilty_dict[string_time]
			run_flag=False
		else:
			for time_key in chunk_feasibilty_dict:
				if abs(float(time_key)-current_time)<1e-5 and media_non_photon_check_flag:
					run_info=chunk_feasibilty_dict[time_key]
					run_flag=False
		print(run_flag)
		if run_flag:
			overall_run_flag=True
			#Adjust the model to the time point
			gene_value_dict=gene_value_dict_maker(float(current_time), total_gene_set, conversion_file, gene_max_dict, model_parameter_file, run_config_dict)
			chunked_model=model_bounds_adjuster(chunked_model,fixed_chunked_model, gene_value_dict, reaction_gene_dict, max_flux, float(current_time), run_config_dict)

			if not media_non_photon_check_flag:
				#If we need to overwrite the media
				print('Media Difference, overwriting {}'.format(string_time))
				if string_time in chunk_feasibilty_dict:
					print(chunk_feasibilty_dict[string_time]['Exchange Bounds']==ex_dict)
					pprint(chunk_feasibilty_dict[string_time]['Exchange Bounds'])
					pprint(ex_dict)
				#overwrite the exchange dict
				chunk_feasibilty_dict[string_time]={}
				print()
				print(string_time)
				chunk_feasibilty_dict[string_time]['Exchange Bounds']=ex_dict

			aa_dict={"ala_c":"ala_L_c",
			"arg_c":"arg_L_c",
			"asn_c":"asn_L_c",
			"asp_c":"asp_L_c",
			"cys_c":"cys_L_c",
			"gln_c":"gln_L_c",
			"glu_c":"glu_L_c",
			"gly_c":"gly_c",
			"his_c":"his_L_c",
			"ile_c":"ile_L_c",
			"leu_c":"leu_L_c",
			"lys_c":"lys_L_c",
			"met_c":"met_L_c",
			"phe_c":"phe_L_c",
			"pro_c":"pro_L_c",
			"ser_c":"ser_L_c",
			"thr_c":"thr_L_c",
			"trp_c":"trp_L_c",
			"tyr_c":"tyr_L_c",
			"val_c":"val_L_c"}
			meta_export_dict={}

			for met in met_list:
				if met not in chunk_feasibilty_dict[string_time]:
					chunk_feasibilty_dict[string_time][met]={}
					with chunked_model as testing_model:
						print('{:25s}'.format(met),end='')
						chunk_required_export_reaction_set=set([])
						goal_met_rxn=testing_model.reactions.get_by_id('{}_formation'.format(met))
						testing_model.objective={goal_met_rxn:1}
						output=testing_model.slim_optimize()
						#print(met)
						#print(output)
						chunk_feasibilty_dict[string_time][met]['Media']=abs(output) if abs(output) > 1e-10 else 0
						print('{}'.format('T' if abs(output) > 1e-10 else 'F'),end=' ')

						for met2 in met_list:
							if met!=met2:
								with testing_model as meta_testing_model:
									meta_testing_model.reactions.get_by_id( '{}_export'.format(met2)).bounds=(-1,1000)
									form_rxn=meta_testing_model.reactions.get_by_id( '{}_formation'.format(met2))
									meta_ex_meta_list=[] #adjust export reactions for metabolites
									for meta_ex_meta in form_rxn.metabolites:
										if 'CH_' not in meta_ex_meta.id:
											meta_factor=-form_rxn.metabolites[meta_ex_meta]
											if 'trna' in meta_ex_meta.id: #trna correction
												if meta_ex_meta.id.startswith('trna'):
													meta_ex_meta=meta_testing_model.metabolites.get_by_id(aa_dict[meta_ex_meta.id[4:]])
													meta_factor*=-1
												else:
													continue

											meta_ex_rxn=meta_testing_model.reactions.get_by_id('{}_chunk_export'.format(meta_ex_meta))
											meta_ex_rxn.bounds=(0, meta_factor)
									print('[',end='',flush=True)
									#with timeout(seconds=3):
									sol=testing_model.optimize()
									output=sol.objective_value
									chunk_feasibilty_dict[string_time][met][met2]=abs(output) if abs(output) > 1e-10 else 0
									#print('{}'.format('T' if abs(output) > 1e-10 else 'F'),end='')
									if abs(output)>1:
										print_out='+'
									elif abs(output)<1e-10:
										print_out='F'
									else:
										print_out=int(math.floor(abs(output)*10))
									print(print_out,end='')
									print(']',end=' ',flush=True)
									#print()
									'''
									if met2 not in meta_export_dict:
										meta_export_dict[met2]=set([])
									for meta_ex_meta in meta_ex_meta_list:
										if abs(sol.shadow_prices[meta_ex_meta.id])>1e-6:
											meta_export_dict[met2].add(meta_ex_meta.id)
									'''

							else:
								print('[S]',end=' ',flush=True)

						print()
						#test shadow prices

						atp_NGAM_rxn=testing_model.reactions.get_by_id('ATPM_NGAM')
						atp_NGAM_rxn.bounds=(atp_NGAM_rxn.bounds[0],1000)
						testing_model.objective={atp_NGAM_rxn:1}
						sol=testing_model.optimize()
						base_ATP=atp_NGAM_rxn.flux
						form_rxn=testing_model.reactions.get_by_id('{}_formation'.format(met))
						chunk_feasibilty_dict[string_time][met]['Export Set']=set([])
						shadow_meta_list=[]
						#print((string_time, met))
						#pprint(base_chunk_feasibility_dict[string_time][met])
						#print(chunk_feasibilty_dict[string_time][met]['Export Set'])
						if string_time in base_chunk_feasibility_dict:
							if met in base_chunk_feasibility_dict[string_time]:
								if 'Export Set' in base_chunk_feasibility_dict[string_time][met]:
									for rxn_id in base_chunk_feasibility_dict[string_time][met]['Export Set']:
											chunk_feasibilty_dict[string_time][met]['Export Set'].add(rxn_id)
						#print(chunk_feasibilty_dict[string_time][met]['Export Set'])
						#sys.exit()
						for meta_ex_meta in form_rxn.metabolites:
							if 'CH_' not in meta_ex_meta.id:
								if 'trna' in meta_ex_meta.id: #trna correction
									if meta_ex_meta.id.startswith('trna'):
										meta_ex_meta=testing_model.metabolites.get_by_id(aa_dict[meta_ex_meta.id[4:]])
									else:
										continue
								shadow_meta_list.append(meta_ex_meta)
								meta_ex_rxn=testing_model.reactions.get_by_id('{}_chunk_export'.format(meta_ex_meta))
								if meta_ex_rxn.id not in chunk_feasibilty_dict[string_time][met]['Export Set']:
									meta_ex_rxn.bounds=(-1, 0)
									#meta_ex_meta_list.append(meta_ex_meta)
									sol=testing_model.optimize()
									new_ATP=atp_NGAM_rxn.flux
									#print(meta_ex_meta, new_ATP, new_ATP-base_ATP)
									if new_ATP-base_ATP<1e-6:
										chunk_feasibilty_dict[string_time][met]['Export Set'].add(meta_ex_rxn.id)
									meta_ex_rxn.bounds=(0, 0)
								#if met in ['CH_Nucleotide']:
								#	print(meta_ex_meta)
								#	print(base_ATP)
								#	print(new_ATP)
								#	print()
								#	pprint(chunk_feasibilty_dict[string_time][met]['Export Set'])


						'''
						shadow_set=set([])
						form_rxn=testing_model.reactions.get_by_id('{}_formation'.format(met))
						ex_rxn=testing_model.reactions.get_by_id('{}_export'.format(met))
						form_rxn.bounds=(-1,1000)
						ex_rxn.bounds=(-1,1000)
						atp_NGAM_rxn.bounds=(atp_NGAM_rxn.bounds[0],1000)
						testing_model.objective={atp_NGAM_rxn:1}
						sol=testing_model.optimize()
						base_ATP=atp_NGAM_rxn.flux
						for meta_ex_meta in shadow_meta_list:
							if abs(sol.shadow_prices[meta_ex_meta.id])<1e-6:
								shadow_set.add('{}_chunk_export'.format(meta_ex_meta))

						pprint(shadow_set)
						print()
						pprint(chunk_feasibilty_dict[string_time][met]['Export Set'])
						print()

						assert shadow_set==chunk_feasibilty_dict[string_time][met]['Export Set']
						'''

						chunk_feasibilty_dict[string_time][met]['Export Set']=list(chunk_feasibilty_dict[string_time][met]['Export Set'])
			#sys.exit()

						#pprint(chunk_feasibilty_dict[string_time][met]['Export Set'])
			#pprint(chunk_feasibilty_dict[string_time])
			#sys.exit()
			#sys.exit()
	return chunk_feasibilty_dict,overall_run_flag

def model_degradation_priority_assessor(run_config_dict, chunk_track_dict, chunk_degradation_error_dict, run_info, chunked_model, mass_dict, debug=False):
	#pass
	#this function decides which chunked metabolites should be degraded
	debug=run_config_dict["Debug"]

	#first assess which chunks can be degraded
	chunk_degradation_set=set([])
	#pprint(run_info)
	#sys.exit()
	for key in run_config_dict['Chunk P Degredation Value']:
		#print(key)
		#print(len(run_info[key]['Export Set'])<(len(chunked_model.reactions.get_by_id('{}_formation'.format(key)).metabolites)-1))
		#print((run_config_dict['Chunk P Degredation Value'][key]!=0 or run_config_dict['Chunk I Degredation Value'][key]!=0))
		if (run_config_dict['Chunk P Degredation Value'][key]!=0 or run_config_dict['Chunk I Degredation Value'][key]!=0) and chunk_track_dict[key]>min(0,run_config_dict['Chunk Bounds Fractions'][key][0]):

			degradation_flag=len(run_info[key]['Export Set'])<(len(chunked_model.reactions.get_by_id('{}_formation'.format(key)).metabolites)-1)

			for check_meta in run_info:
				if 'CH' in check_meta:
					if run_info[check_meta]['Media']:
						degradation_flag=True
						break
			if not degradation_flag:
				for check_meta in run_info:
					if 'CH' in check_meta:
						if key in run_info[check_meta]:
							if run_info[check_meta][key]:
								degradation_flag=True
								break
			if degradation_flag:
				chunk_degradation_set.add(key)

	#next, find the total error. Error is positive when the chunk is above the setpoint (though integral error terms can adjust this).
	total_degradation_error=0
	for key in chunk_degradation_set:
		total_degradation_error+=chunk_degradation_error_dict[key]

	if debug:
		print('total_degradation_error')
		print(total_degradation_error)
		print()

	#finally, calculate the weight and sort it
	deg_chunk_list=[]
	deg_value_list=[]
	for key in chunk_degradation_set:
		mass_diff=0
		for reaction_id in run_info[key]['Export Set']:
			meta_form_reaction=chunked_model.reactions.get_by_id('{}_formation'.format(key))
			if chunked_model.metabolites.get_by_id(reaction_id[:-13]) in meta_form_reaction.metabolites:
				factor=-1 * meta_form_reaction.metabolites[chunked_model.metabolites.get_by_id(reaction_id[:-13])]
			else: #trna correction
				factor=meta_form_reaction.metabolites[chunked_model.metabolites.get_by_id('trna{}_c'.format(reaction_id[:3]))]
			for ex_meta in chunked_model.reactions.get_by_id(reaction_id).metabolites:
				meta_factor=-chunked_model.reactions.get_by_id(reaction_id).metabolites[ex_meta]
			mass_diff+=factor/meta_factor
		deg_chunk_list.append(key)
		deg_value_list.append(mass_diff*chunk_degradation_error_dict[key]/abs(total_degradation_error))

	deg_keydict = dict(zip(deg_chunk_list, deg_value_list))
	deg_chunk_list.sort(key=deg_keydict.get,reverse=True)
	if debug:
		print(deg_chunk_list, deg_keydict)
		sys.exit()
	return deg_chunk_list, deg_keydict


def model_production_priority_assessor(run_config_dict, chunk_track_dict, chunk_production_error_dict, run_info, deg_chunk_list,debug=False):
	#pass
	#this function decides which chunked metabolites should be produced
	#NOTE it only assesses chunks with positive error

	debug=run_config_dict["Debug"]

	#first, assess which metabolites can be produced
	chunk_production_set=set([])
	for key in run_info:
		if 'CH' in key:
			for meta in run_info[key]:
				if run_info[key][meta] and (meta in deg_chunk_list or meta=='Media'):
					chunk_production_set.add(key)
					break
	if debug:
		print('chunk_production_set')
		print(chunk_production_set)
		print()
	#second, calculate the weight of the error of the metabolites that can be produced
	total_production_error=0
	for key in chunk_production_set:
		if chunk_production_error_dict[key]>0: #a negative production error means that the metabolite in question is above the setpoint
			total_production_error+=chunk_production_error_dict[key]
	if debug:
		print('total_production_error')
		print(total_production_error)
		print()
	if total_production_error==0: #if all of these are under the error, produce all of them at the ratio defined by the chunk goal
		total_production_error=1
		for key in run_config_dict['Chunk Goal Information']:
			chunk_production_error_dict[key]=run_config_dict['Chunk Goal Information'][key]
	if debug:
		print('chunk_production_error_dict')
		pprint(chunk_production_error_dict)
		print()
	#third, calculate the weight and sort
	prod_chunk_list=[]
	prod_value_list=[]
	for key in chunk_production_set:
		if chunk_production_error_dict[key]>0:
			prod_chunk_list.append(key)
			prod_value_list.append(chunk_production_error_dict[key]/total_production_error)

	prod_keydict = dict(zip(prod_chunk_list, prod_value_list))
	prod_chunk_list.sort(key=prod_keydict.get, reverse=True)
	return prod_chunk_list, prod_keydict



def shadow_production_priority_assessor(run_config_dict, chunk_track_dict, chunk_production_error_dict, run_info, deg_chunk_list,debug=False):
	#pass
	#this function decides which chunked metabolites should be produced
	#NOTE it only assesses chunks with positive error

	debug=run_config_dict["Debug"]

	#first, assess which metabolites can be produced
	chunk_production_set=set([])
	total_mass=0
	for meta in chunk_track_dict:
		total_mass+=chunk_track_dict[meta]
	for key in run_info:
		if 'CH' in key and key!='CH_Energy' and chunk_track_dict[key]/total_mass < run_config_dict['Chunk Bounds Fractions'][key][1]:
			for meta in run_info[key]:
				if run_info[key][meta] and (meta in deg_chunk_list or meta=='Media'):
					chunk_production_set.add(key)
					break
	if debug:
		print('chunk_production_set')
		print(chunk_production_set)
		print()

	#pprint(chunk_production_error_dict)
	#map error to exponential space
	total_err_val=0
	prod_chunk_list=[]
	prod_keydict={}
	for key in chunk_production_set:
		prod_chunk_list.append(key)
		prod_keydict[key]=10**chunk_production_error_dict[key] #to map the distances into weighting 0->inf space
		total_err_val+=prod_keydict[key]

	for key in prod_keydict:
		prod_keydict[key]/=total_err_val #normalize to 1

	prod_chunk_list.sort(key=prod_keydict.get, reverse=True)
	return prod_chunk_list, prod_keydict




def model_time_evaluator(run_config_dict, chunk_track_dict, chunk_production_error_dict, chunk_degradation_error_dict, chunked_model, fixed_chunked_model, total_gene_set, conversion_file, gene_max_dict, reaction_gene_dict, chunk_feasibilty_dict, sim_time, time_interval, flux_export_dict, mass_ex_min_obj_obj, mass_dict):
	#This function evaluates the model at a specific time
	#first, check the time to see if the feasibility already exists
	#print('\n\n\n')
	current_time=sim_time%24
	run_flag=True
	string_time=time_stringer(current_time)

	if string_time in chunk_feasibilty_dict:
		run_info=chunk_feasibilty_dict[string_time]
		run_flag=False
	else:
		for current_time_key in chunk_feasibilty_dict:
			if abs(float(current_time_key)-current_time)<1e-4:
				run_info=chunk_feasibilty_dict[current_time_key]
				run_flag=False
	if run_flag:
		print('Time data not available, running {}'.format(current_time))
		chunk_feasibilty_dict, overall_run_flag=model_feasibility_tester(chunked_model,fixed_chunked_model, total_gene_set, conversion_file, gene_max_dict, reaction_gene_dict, chunk_feasibilty_dict, run_config_dict, [current_time])
		run_info=chunk_feasibilty_dict[string_time]

	print('C',end='',flush=True)
	deg_chunk_list, deg_keydict=model_degradation_priority_assessor(run_config_dict, chunk_track_dict, chunk_degradation_error_dict, run_info, chunked_model, mass_dict)
	print('D',end='',flush=True)
	prod_chunk_list, prod_keydict=shadow_production_priority_assessor(run_config_dict, chunk_track_dict, chunk_production_error_dict, run_info, deg_chunk_list)
	print('P',end='',flush=True)
	#print()
	#pprint(deg_keydict)
	#pprint(deg_chunk_list)
	#print()

	#print()
	#pprint(prod_keydict)
	#pprint(prod_chunk_list)
	#print()

	gene_value_dict=gene_value_dict_maker(float(current_time), total_gene_set, conversion_file, gene_max_dict,model_parameter_file, run_config_dict)
	chunked_model=model_bounds_adjuster(chunked_model,fixed_chunked_model, gene_value_dict, reaction_gene_dict, max_flux, float(current_time), run_config_dict)


	#first, make NGAM
	with chunked_model as model:
		light_rxn=model.reactions.get_by_id('EX_photonVis_e')
		print('LLL:{0:6.1f} '.format(light_rxn.bounds[0]),end='')
		atp_NGAM_rxn=model.reactions.get_by_id('ATPM_NGAM')
		ngam_bounds=atp_NGAM_rxn.bounds
		upper_NGAM=ngam_bounds[1]
		atp_NGAM_rxn.bounds=[0,ngam_bounds[1]]
		model.objective={atp_NGAM_rxn:1}
		ngam_output=model.slim_optimize()
		mediaFlag= (not np.isnan(ngam_output)) and ngam_output>=ngam_bounds[1]
		#if NGAM can't be made with media, check the consumption
		if mediaFlag:
			print('M ',end='\t')
		else:
			print('!M',end='\t')

		#print(model.reactions.CH_Energy_export.bounds)
		model.reactions.DM_na1_h.bounds=(0,0)
		#print(model.reactions.CH_Energy_export.bounds)

		#for reaction in model.reactions:
		#	if 'chunk_export' in reaction.id:
		#		print(reaction)
		#		reaction.bounds=(0,1000)


		ngam_lower_bound=0
		index=0
		check_val_list=[]
		if bool(run_config_dict['Debug']):
			print('\ndeg_chunk_list')
			pprint(deg_chunk_list)
			print()
		#While ngam is not reaching the max, check the possible catabolic sources
		#print(deg_chunk_list)

		if (np.isnan(ngam_output) or ngam_output<ngam_bounds[1]) and len(deg_chunk_list)>0:
			chunk_mass_constraint=None
			for meta in deg_chunk_list:
				meta_ex_reaction=model.reactions.get_by_id('{}_export'.format(meta))
				#print(meta_ex_reaction)
				#print(meta_ex_reaction.flux_expression)
				if chunk_mass_constraint==None:
					chunk_mass_constraint=-1*meta_ex_reaction.reverse_variable
				else:
					chunk_mass_constraint+=-1*meta_ex_reaction.reverse_variable
		#print()
		#print(chunk_mass_constraint)
		#print()
			mass_ex=model.problem.Constraint(chunk_mass_constraint, lb=-run_config_dict['Max Degradation Rate'], ub=run_config_dict['Max_flux'], name='Chunk_Degradation_Cons')
		#print(mass_ex)
			consname=model.add_cons_vars(mass_ex)
			print('𝛘',end='',flush=True)

		chunk_mass_min_objective=None
		meta_rxn_list=[]
		total_mass=0
		for key in chunk_track_dict:
			total_mass+=chunk_track_dict[key]

		while (np.isnan(ngam_output) or ngam_output<ngam_bounds[1]) and index<len(deg_chunk_list):
			#go in order, starting with the first (which has the best score, essentially)
			meta=deg_chunk_list[index]
			#print(meta)
			meta_ex_reaction=model.reactions.get_by_id('{}_export'.format(meta))
			meta_form_reaction=model.reactions.get_by_id('{}_formation'.format(meta))
			#print(meta_ex_reaction)
			#print(meta_ex_reaction.bounds)
			meta_ex_reaction.bounds=(-min(run_config_dict['Max Degradation Rate'], chunk_track_dict[meta]/total_mass), meta_ex_reaction.bounds[1]) #Fix the exchange reaction to the max or whatever remains, and the previous dump rate
			#print()
			#print(meta)
			#print(run_info[meta]['Export Set'])
			#print()
			mass_diff=0
			for reaction_id in run_info[meta]['Export Set']:
				if model.metabolites.get_by_id(reaction_id[:-13]) in meta_form_reaction.metabolites:
					factor=-1 * meta_form_reaction.metabolites[model.metabolites.get_by_id(reaction_id[:-13])]
				else: #trna correction
					factor=meta_form_reaction.metabolites[model.metabolites.get_by_id('trna{}_c'.format(reaction_id[:3]))]
				for ex_meta in model.reactions.get_by_id(reaction_id).metabolites:
					meta_factor=-model.reactions.get_by_id(reaction_id).metabolites[ex_meta]
				mass_diff+=factor/meta_factor
				#test

				factor=1
				meta_factor=1
				model.reactions.get_by_id(reaction_id).bounds=(0, run_config_dict['Max Degradation Rate']*factor/meta_factor)

				#print(model.reactions.get_by_id(reaction_id), model.reactions.get_by_id(reaction_id).bounds)
			#print(meta+' Mass Immediately Retained', mass_diff)
			if chunk_mass_min_objective==None:
				chunk_mass_min_objective=mass_diff*meta_ex_reaction.reverse_variable
			else:
				chunk_mass_min_objective+=mass_diff*meta_ex_reaction.reverse_variable

			#print(meta_ex_reaction)
			#print(meta_ex_reaction.bounds)
			#psol=cobra.flux_analysis.pfba(model)
			model.slim_optimize()
			ngam_output=atp_NGAM_rxn.flux
			#if meta=='CH_Nucleotide':
			#	print(atp_NGAM_rxn.flux, atp_NGAM_rxn.bounds, atp_NGAM_rxn)
			#	print(meta_ex_reaction.flux, meta_ex_reaction.bounds, meta_ex_reaction)
			#print((meta_ex_reaction.flux, ngam_output))
			# if there is an improvement in the ngam production, adjust it.
			if not np.isnan(ngam_output) and ngam_output>ngam_lower_bound:
				ngam_lower_bound=min(ngam_output, ngam_bounds[1])
				atp_NGAM_rxn.bounds=[ngam_lower_bound,ngam_bounds[1]]
				#meta_ex_reaction.bounds=(-run_config_dict['Max Degradation Rate'], max(meta_ex_reaction.flux, -run_config_dict['Max Degradation Rate']))


			#print(meta_ex_reaction.bounds)
			check_val_list.append((meta,meta_ex_reaction.id))
			#for meta, rxn_id in check_val_list:
			#	print(meta, model.reactions.get_by_id(rxn_id).bounds)
			index+=1
			#print()
			#for reaction in model.reactions:
			#	if len(reaction.metabolites)==1 and abs(reaction.flux)>1e-6:
			#		print(reaction, reaction.flux)
			#print('\n'*4)
			meta_rxn_list.append(meta_ex_reaction)

		#print(chunk_mass_min_objective)
		if atp_NGAM_rxn.flux<upper_NGAM and False:
			for reaction in model.reactions:
				if 'CH_' in reaction.id and 'export' in reaction.id and abs(reaction.flux)>1e-6:
					print(reaction.flux, reaction.bounds, reaction)
			print()
			for reaction in model.reactions:
				if 'CH_' in reaction.id and 'formation' in reaction.id and abs(reaction.flux)>1e-6:
					print(reaction.flux, reaction.bounds, reaction)
			print()
			for reaction in model.reactions:
				if abs(reaction.flux)>1e-8 and (abs(reaction.flux-reaction.bounds[0])<1e-8 or abs(reaction.flux-reaction.bounds[1])<1e-8):
					print(reaction.flux, reaction.bounds, reaction)
			sys.exit()


		#mass_ex=model.problem.Constraint(chunk_mass_constraint, lb=-run_config_dict['Max Degradation Rate'], ub=run_config_dict['Max_flux'], name='Chunk_Degradation_Cons')

		if bool(run_config_dict['Debug']) or False:
			print(model.constraints['Chunk_Degradation_Cons'])
			for rxn in meta_rxn_list:
				print(rxn)
				print(rxn.bounds)
				print(rxn.forward_variable.primal)
				print(rxn.reverse_variable.primal)
				print()

			print(atp_NGAM_rxn.flux)
			print(atp_NGAM_rxn.bounds)
			#sys.exit()
			print('Exchange Reaction Breakdown')
			for reaction in model.reactions:
				if '_chunk_export' in reaction.id and False:
					print(reaction.flux, reaction.bounds, reaction)
					if abs(reaction.flux)>1e-6:
						for metabolite in reaction.metabolites:
							for meta_reaction in metabolite.reactions:
								if abs(meta_reaction.flux)>1e-6:
									print('\t',end='')
									print(meta_reaction.flux, meta_reaction.bounds, meta_reaction)

		#sys.exit()


		#sys.exit()
		print('NGAM:{:4.3f}'.format(atp_NGAM_rxn.flux),end='\t',flush=True)
		#print(model.reactions.CH_Energy_export.bounds)
		#print(model.reactions.CH_Energy_export.flux)


		if False:
			print()
			print('\n'*10)
			pprint(model.medium)
			print(model.reactions.EX_photonVis_e)
			print(model.reactions.EX_photonVis_e.flux)
			print(atp_NGAM_rxn)
			print(atp_NGAM_rxn.flux)
			for exchange in model.exchanges:
				if abs(exchange.flux)>1e-6:
					print(exchange)
					print(exchange.flux)
					print(exchange.bounds)
			#sys.exit()

			print(model.reactions.EX_photonVis_e)
			print(model.reactions.EX_photonVis_e.flux)
			for exchange in model.exchanges:
				if abs(exchange.flux)>1e-6:
					print(exchange)
					print(exchange.flux)
					print(exchange.bounds)
			print(prod_chunk_list)

			for meta in prod_chunk_list:
				goal_rxn=model.reactions.get_by_id('{}_export'.format(meta))
				print(goal_rxn)
				print(goal_rxn.bounds)
				model.objective={goal_rxn:1}
				model.optimize()
				print(goal_rxn.flux)

			#sys.exit()


		for exchange in model.exchanges:
			if exchange.flux<exchange.bounds[0]-1e-6 or exchange.flux>exchange.bounds[1]+1e-6:
				print('ERROR')
				print(exchange)
				print(exchange.flux)
				print(exchange.bounds)
				sys.exit()
		if bool(run_config_dict['Debug']):
			for reaction in model.reactions:
				if reaction.flux<reaction.bounds[0]-1e-6 or reaction.flux>reaction.bounds[1]+1e-6:
					print('ERROR')
					print(reaction)
					print(reaction.flux)
					print(reaction.bounds)
					sys.exit()
		#Next, check the production status. If all metabolites with a positive error term have been drawn down, then don't produce more.
		index-=1 #adjust for the end of the while loop increment
		produceFlag=True
		if not mediaFlag:
			produceFlag=False
			if ngam_output<ngam_bounds[1]: #if the cell cannot make maintenance ATP it is obviously in survival mode
				produceFlag=False
			for rxn in meta_rxn_list: #If the cell is drawing from negative stores to stay alive
				if rxn.bounds[0]<-1e-6:
					if rxn.id[:-7] not in deg_keydict:
						print(rxn)
						pprint(deg_keydict)
						print(deg_chunk_list)
					if deg_keydict[rxn.id[:-7]]<0:
						produceFlag=False
		#print(model.reactions.CH_Energy_export.bounds)


		print('PR:{}'.format(int(produceFlag)),end=' ', flush=True)
		#sys.exit()

		meta_dict={'CH_Nucleotide_c':'Nu', 'CH_Lipid_c':'Li', 'CH_LipidDroplet_c':'LD', 'CH_Carbohydrate_c':'Crb', 'CH_Carotenoid_c':'Cro', 'CH_Protein_c':'Pt','CH_Starch_c':'St','CH_Redox_c':'Rdx', 'CH_Chla_c':'Ca', 'CH_Chlb_c':'Cb'}
		with model as product_model:
			product_flux=np.inf
			if produceFlag:
				prev_prod_list=None
				#print('Here')
				prod_index=0
				#knockout exchanges
				for meta in prod_chunk_list:
					meta_bounds=product_model.reactions.get_by_id('{}_export'.format(meta)).bounds
					new_lb=meta_bounds[0]
					if meta_bounds[1]<0:
						new_ub=meta_bounds[1]
					else:
						new_ub=0
					#product_model.reactions.get_by_id('{}_export'.format(meta)).bounds= (new_lb,new_ub)
				#product_model.reactions.CH_Energy_export.bounds=(0,1000)
				#repeatedly check the reactions until the shadow prices of each chunk is non-zero
				prod_list=prod_chunk_list
				iter_count=0
				production_reaction_dict={}
				prev_dict=None
				prev2_dict=None
				while len(prod_list)>0 and prev_prod_list!=prod_list and production_reaction_dict!=prev2_dict:

					#print(prod_list)
					prod_coeff_sum=0
					for meta in prod_list:
						prod_coeff_sum+=prod_keydict[meta]
					#assert abs(prod_coeff_sum-1)<1e-4

					#create a new production reaction
					prev2_dict=prev_dict
					prev_dict=production_reaction_dict
					production_reaction=cobra.Reaction('product_reaction_{}'.format(iter_count))
					form_ATP=run_config_dict['Growth ATP']
					production_reaction_dict={'atp_c':-form_ATP,'h2o_c':-form_ATP,'adp_c':form_ATP,'h_c':form_ATP,'pi_c':form_ATP}
					new_prod_coeff_sum=0
					final_prod_list=[]
					for meta in prod_list:
						#remove small values (speeds up optimization and objective finding)
						if prod_keydict[meta]/prod_coeff_sum>1e-2:
							new_prod_coeff_sum+=prod_keydict[meta]
							final_prod_list.append(meta)
					for meta in final_prod_list:
						production_weight=prod_keydict[meta]/new_prod_coeff_sum #normalize so that all fluxes that go through production reaction have a 1 to 1 ratio with chunk
						#if True: #prevents ping-pong behavior
							#production_weight=production_weight**2/(production_weight+1e-4)
						#	production_weight=np.sqrt(max(production_weight**2-0.01,0))
						production_reaction_dict['{}_c'.format(meta)]=-production_weight
					product_model.add_reaction(production_reaction)
					production_reaction.bounds=(0,1000)
					#print('addedrxn')
					if False:
						print()
						pprint(production_reaction_dict)
						print()
					prev_dict=production_reaction_dict
					production_reaction.add_metabolites(production_reaction_dict)
					#print(production_reaction)
					product_model.objective=production_reaction
					#print('{0}'.format(iter_count),end='|')
					if False:
						print()
						print()

						for meta in prod_list:
							meta_formation_rxn= product_model.reactions.get_by_id('{}_formation'.format(meta))
							meta_export_rxn= product_model.reactions.get_by_id('{}_export'.format(meta))
							product_model.objective={meta_export_rxn:1}
							print(meta_export_rxn)
							print(meta_export_rxn.bounds)
							print(meta_formation_rxn.id)
							print(meta_formation_rxn.bounds)
							product_model.optimize()
							print(meta_export_rxn.flux)
							print(meta_formation_rxn.flux)

						export_test=cobra.Reaction('TEST_Export')
						export_test.name='TEST OF EXPORTING'
						model.add_reactions([export_test])
						export_test.add_metabolites({'atp_c':-1,'h2o_c':-1,'adp_c':1,'h_c':1,'pi_c':1})
						model.objective={export_test:1}
						model.optimize()
						print()

						for reaction in model.reactions:
							if abs(reaction.flux)>1e-6:
								print(reaction)
								print(reaction.flux)
								print(reaction.bounds)

						print()
						print(export_test)
						print(export_test.flux)

						print('\n\n')
						pprint(production_reaction_dict)
						print(production_reaction)
						print(production_reaction.bounds)
						product_model.objective={production_reaction:1}
						product_model.optimize()
						print(production_reaction.flux)
						sys.exit()

					#maximize flux through this production reaction
					solution=product_model.optimize()
					#solution=cobra.flux_analysis.pfba(product_model)
					product_flux=production_reaction.flux
					if np.isnan(product_flux):
						print('here')
						product_flux=0
					#fix this production
					production_reaction.bounds=(product_flux,product_flux)
					if abs(product_flux)>1e-8:
						if iter_count>0:
							print('',flush=True)
						print('{0:3.2}'.format(prod_coeff_sum),end=':')
						print('{0:5.2}'.format(product_flux),end=', ', flush=True)
						for meta in production_reaction.metabolites:
							if 'CH_' in meta.id:
								print('{0}:{1:5.3f};'.format(meta_dict[meta.id], -production_reaction.metabolites[meta]),end='')
						print(', ',end='', flush=True)
						#print()
						#print(production_reaction)
						#sys.exit()
						if False:
							print()
							for meta in production_reaction_dict:
								print('{} : {}'.format(meta, production_reaction_dict[meta]))
							print()
							print(solution.shadow_prices)
							for metabolite in production_reaction.metabolites:
								print('{} : {}'.format(metabolite, solution.shadow_prices[metabolite.id]))
							print('{} : {}'.format('nh4_e', solution.shadow_prices['nh4_e']))
							#sys.exit()
							for exchange in product_model.exchanges:
								if abs(exchange.flux)>1e-6:
									print(exchange)
									print(exchange.flux)
									print()

							for reaction in product_model.reactions:
								if ('chunk_export' in reaction.id or 'CH_' in reaction.id) and abs(reaction.flux)>1e-7:
									print(reaction.id)
									print(reaction)
									print(reaction.flux)
									print()
							print()
							#sys.exit()

					#evaluate shadow prices
					prev_prod_list=final_prod_list
					for meta in final_prod_list:
						#print(meta, solution.shadow_prices['{}_c'.format(meta)])
						if solution.shadow_prices['{}_c'.format(meta)]<-1e-8:
							prod_list.remove(meta)
					#print(prod_list)
					iter_count+=1



				#pprint(chunk_track_dict)
			elif abs(atp_NGAM_rxn.flux)>1e-3:
				ch_ex_obj=model.problem.Objective(chunk_mass_min_objective, direction='min')
				model.objective=ch_ex_obj
				sol=model.optimize(objective_sense=None)
				try:
					sol=cobra.flux_analysis.pfba(model)
				except:
					pass
			print('',end='\t',flush=True)
			#update model values for valid cells
			if not np.isnan(ngam_output):
				#Print the usage rates of each consumed metabolite
				print('C_V:',end=' ')
				for item in check_val_list:
					print('{0}:{1:4.4f};'.format(meta_dict[item[0]+'_c'], model.reactions.get_by_id(item[1]).flux),end='')
				print('',end='\t')

				failed_run_Flag=False
				for reaction in product_model.reactions:
					if reaction.flux<reaction.bounds[0]-1e-4 or reaction.flux>reaction.bounds[1]+1e-4:
						print(reaction)
						print(reaction.flux, reaction.bounds)
						failed_run_Flag=True

				if failed_run_Flag:
					sys.exit()

				if False and (sim_time>12 and abs(product_model.reactions.CH_Carbohydrate_export.flux)>1e-2):
					cobra.io.write_sbml_model(product_model,'Time_Specific_Model.xml')
					print('\n\n')
					for reaction in model.reactions:
						if len(reaction.metabolites)==1 and abs(reaction.flux)>1e-4:
							print(reaction.flux, reaction.bounds, reaction)
					print('\n\n')
					for reaction in model.reactions:
						if abs(reaction.flux)>1e-4 and (abs(reaction.flux-reaction.bounds[0])<1e-4 or abs(reaction.flux-reaction.bounds[1])<1e-4):
							print(reaction.flux, reaction.bounds, reaction)

					sys.exit()
				'''
				try:
					product_model.objective=mass_ex_min_obj_obj
					product_model.objective={atp_NGAM_rxn:1}
					cobra.flux_analysis.pfba(product_model)

				except:
					print('Min mass ex failed')
					for reaction in product_model.reactions:
						if reaction.flux<reaction.bounds[0] or reaction.flux>reaction.bounds[1]:
							print(reaction)
							print(reaction.flux, reaction.bounds)
					print()
					sys.exit()
					product_model.objective={atp_NGAM_rxn:1}
					cobra.flux_analysis.pfba(product_model)
				#product_model.optimize()
				'''

				#print('OPTIMIZING')
				#sys.exit()
				#print(product_model.reactions.CH_Energy_export.bounds)
				#print(product_model.reactions.CH_Energy_export.flux)
				total_mass=0
				for meta in chunk_track_dict:
					total_mass+=chunk_track_dict[meta]
				total_exchange_flux=0
				#print()
				#print(atp_NGAM_rxn.flux)
				if False:
					print()
					for exchange in product_model.exchanges:
						if abs(exchange.flux)>1e-6:
							print(exchange)
							print(exchange.flux)
							print(exchange.bounds)
							print()

					for reaction in product_model.reactions:
						if ('chunk_export' in reaction.id or 'CH_' in reaction.id) and abs(reaction.flux)>1e-7:
							print(reaction.id)
							print(reaction)
							print(reaction.flux)
							print(reaction.bounds)
							print()
					print()
				for meta in chunk_track_dict:
					if meta != 'CH_Energy':
						meta_reaction=product_model.reactions.get_by_id('{}_formation'.format(meta))
						#if check_val<0:
							#print((meta, meta_reaction.flux*time_interval))
						chunk_track_dict[meta]+=meta_reaction.flux*time_interval*total_mass
						chunk_track_dict[meta]=max(chunk_track_dict[meta],0)
						total_exchange_flux+=meta_reaction.flux*time_interval*total_mass
					#pprint(chunk_track_dict)
				print('TEF:{0:6.4f}'.format(total_exchange_flux),end='\t')

				for exchange in model.exchanges:
					if exchange.flux<exchange.bounds[0]-1e-6 or exchange.flux>exchange.bounds[1]+1e-6:
						print('ERROR')
						print(exchange)
						print(exchange.flux)
						print(exchange.bounds)
						sys.exit()

				if (abs(total_exchange_flux)<0.0001 and produceFlag and not mediaFlag) or False:
					print('Error here')
					for exchange in product_model.exchanges:
						if abs(exchange.flux)>1e-6:
							print(exchange)
							print(exchange.bounds)
							print(exchange.flux)
							if exchange.flux<exchange.bounds[0]:
								print('ERROR')
								print('"{}":1e-8'.format(exchange.id))
					cobra.io.save_json_model(product_model,'product_model.json')
					flux_dump_dict={}
					for reaction in product_model.reactions:
						flux_dump_dict[reaction.id]=reaction.flux
					json.dump(flux_dump_dict,open('flux_dump.json','w'))
					cobra.io.write_sbml_model(product_model,'product_model.xml')
					sys.exit()


				light_rxn=product_model.reactions.get_by_id('EX_photonVis_e')
				print('LU: {0:6.1f}'.format(light_rxn.flux),end='\t')

				if abs(total_exchange_flux)<0.0001 and False:
					print()
					print(current_time)
					print(prod_keydict)
					pprint(run_info)
					print()
					print(model_production_priority_assessor(run_config_dict, chunk_track_dict, chunk_production_error_dict, run_info, deg_chunk_list,True))
					sys.exit()


				if flux_export_dict!=None:
					flux_working_dict={}
					bounds_working_dict={}
					for reaction in chunked_model.reactions:
						flux_working_dict[reaction.id]=reaction.flux
						bounds_working_dict[reaction.id]=reaction.bounds
					flux_export_dict['Fluxes'].append(flux_working_dict)
					flux_export_dict['Bounds'].append(bounds_working_dict)
					flux_export_dict['Times'].append(time_stringer(sim_time))

				if run_config_dict['FVA']:
					#print(obligatory_reaction_set)
					obligatory_reactions=json.load(open('{}/FVA/obligatory_reactions.json'.format(storage_directory),'r'))
					print('FVA:',end='',flush=True)
					if time_stringer(sim_time) not in obligatory_reactions:
						FVA_df=flux_variability_analysis(product_model)
						obligatory_reactions[time_stringer(sim_time)]=[]
					#print(FVA_df)
					#print(type(FVA_df))
						for reaction in product_model.reactions:
						#print(reaction)
						#print(FVA_df.loc[reaction.id])
						#print(FVA_df.loc[:,reaction.id])
						#print(FVA_df.loc[reaction.id,:])
						#print(0>FVA_df.loc[reaction.id,'maximum'])
						#print(0<FVA_df.loc[reaction.id,'minimum'])
						#print(0>FVA_df.loc[reaction.id,'maximum'] or 0<FVA_df.loc[reaction.id,'minimum'])
							if (0>FVA_df.loc[reaction.id,'maximum'] or 0<FVA_df.loc[reaction.id,'minimum']):
								obligatory_reactions[time_stringer(sim_time)].append(reaction.id)
						json.dump(obligatory_reactions,open('{}/FVA/obligatory_reactions.json'.format(storage_directory),'w'),indent=4)

					print('C',end='\t',flush=True)
						#sys.exit()
				#pprint(ex_rxn_dict)
				#sys.exit()
				ex_rxn_dict={}
				for ex_rxn in product_model.reactions:
					if len(ex_rxn.metabolites)==1:
						ex_rxn_dict[ex_rxn.id]=ex_rxn.flux

			else:
				flux_export_dict['Fluxes'].append(np.nan)
				bounds_working_dict={}
				for reaction in chunked_model.reactions:
					bounds_working_dict[reaction.id]=reaction.bounds
				flux_export_dict['Bounds'].append(bounds_working_dict)
				#flux_export_dict['Bounds'].append(np.nan)
				ex_rxn_dict={}



	return chunk_feasibilty_dict, chunk_track_dict, flux_export_dict, ex_rxn_dict

def error_calculator(chunk_track_dict, chunk_error_integral_dict, run_config_dict, setpoint_dict):
	#this calculates all the error terms
	chunk_ratio_dict={}
	chunk_ratio_error_dict={}
	#chunk_error_integral_dict={}
	chunk_production_error_dict={}
	chunk_degradation_error_dict={}
	total_chunk_amount=0
	for chunk in chunk_track_dict:
		if chunk != 'CH_Energy':
			total_chunk_amount+=chunk_track_dict[chunk]
	for chunk in chunk_track_dict:
		#chunk_ratio_dict[chunk]=chunk_track_dict[chunk]/total_chunk_amount
		#calculate normalized error. Ignore energy
		if chunk =='CH_Energy':
			chunk_ratio_dict[chunk]=0
			chunk_ratio_error_dict[chunk]=0
		else:
			chunk_ratio_dict[chunk]=chunk_track_dict[chunk]/total_chunk_amount
			chunk_ratio_error_dict[chunk]=(setpoint_dict[chunk]-chunk_ratio_dict[chunk])/setpoint_dict[chunk]
		chunk_error_integral_dict[chunk]+=chunk_ratio_error_dict[chunk]
		chunk_production_error_dict[chunk]=run_config_dict['Chunk P Production Value'][chunk]*chunk_ratio_error_dict[chunk] + run_config_dict['Chunk I Production Value'][chunk]*chunk_error_integral_dict[chunk]
		chunk_degradation_error_dict[chunk]=-(run_config_dict['Chunk P Degredation Value'][chunk]*chunk_ratio_error_dict[chunk] + run_config_dict['Chunk I Degredation Value'][chunk]*chunk_error_integral_dict[chunk])
	#print('Error')
	#pprint(chunk_ratio_error_dict)
	#print()
	return chunk_error_integral_dict, chunk_ratio_dict, chunk_production_error_dict, chunk_degradation_error_dict

def light_and_media_adjuster(run_config_dict, chunked_model, current_time):
	#adjust media as needed
	#pprint(chunked_model.medium)
	input_model_media_dict={}
	for key in run_config_dict['Media']:
		if key[0] is not '#' and key is not 'EX_photonVis_e' and not isinstance(run_config_dict['Media'][key],str):
			input_model_media_dict[key]=run_config_dict['Media'][key]
	chunked_model.medium=input_model_media_dict

	light_info=run_config_dict['Light Regiment']
	light_rxn=chunked_model.reactions.get_by_id('EX_photonVis_e')
	if current_time%light_info['Day']<light_info['Light_Start'] or current_time%light_info['Day']>light_info['Light_End']:
		light_level=0
	else:
		light_level=light_info['Day_Level']
	light_rxn.bounds=(-light_level,light_rxn.bounds[1])
	#print()
	#pprint(chunked_model.medium)
	#sys.exit()

	return chunked_model, light_level

def light_index(storage_dict,light_list):
	index=0
	past_light_level=None
	index_set_list=[]
	r_index=0
	l_index=0
	while index<len(storage_dict['time'])-1:
		if light_list[index]!=past_light_level:
			index_set_list.append([r_index,l_index])
			r_index=l_index
		l_index=index
		index+=1
	return index_set_list

def mass_and_error_plotter(run_config_dict, storage_dict, total_mass_list, light_list, prod_error_storage_dict, deg_error_storage_dict, storage_directory, exchange_storage_dict, setpoint_storage_dict):

	#First, generate the list of light sections
	index_set_list=light_index(storage_dict,light_list)

	fig, ax1=plt.subplots(figsize=(12,8))
	plt.xlabel('Time (Hrs)')
	plt.title('Chunk Mass Ratio')
	for index_set in index_set_list:
		plt.axvspan(storage_dict['time'][index_set[0]],storage_dict['time'][index_set[1]], color='gray', alpha=1-light_list[index_set[0]]/run_config_dict['Light Regiment']['Day_Level'])
	ax1.set_ylabel('Fraction of Biomass')
	for meta in chunk_track_dict:
		if meta in run_config_dict['Chunk Plot']:
			col=run_config_dict['Chunk Color'][meta]
			plt.plot(storage_dict['time'],storage_dict[meta],label=meta, color=col)
			#plt.axhline(y=run_config_dict['Chunk Goal Information'][meta], linestyle='--', xmin=0.85, xmax=1, color=col)
	plt.legend(loc='upper left')
	box = ax1.get_position()
	ax1.set_position([box.x0, box.y0, box.width * 0.7, box.height])

	# Put a legend to the right of the current axis
	ax1.legend(loc='upper left', bbox_to_anchor=(1.1, 0.4))
	ax2=ax1.twinx()
	col='xkcd:black'
	ax2.set_ylabel('Overall Mass')
	ax2.plot(storage_dict['time'],total_mass_list,color=col,label='Total Mass')
	#fig.tight_layout()
	ax2.set_position([box.x0, box.y0, box.width * 0.7, box.height])
	ax2.legend(loc='lower left', bbox_to_anchor=(1.1, 0.8))

	print("Day/Night Masses")
	prev_night=None
	prev_day=None
	for m_time,mass in zip(storage_dict['time'],total_mass_list):
		#print((m_time, abs(m_time%12), mass))
		if m_time%12<0.01 or (12-m_time%12)<0.01:
			if m_time%24<0.01 or (24-m_time%24)<0.01:
				print("Time: {0:8.3f} hrs; mass {1:9.3f}".format(m_time,mass),end='\t')
				if prev_day!=None:
					print("24 hr Ratio: {0:6.3f}".format(mass/prev_day),end='\t')
				prev_day=mass
				print()
			else:
				print("Time: {0:8.3f} hrs; mass {1:9.3f}".format(m_time,mass),end='\t')
				if prev_night!=None:
					print("24 hr Ratio: {0:6.3f}".format(mass/prev_night),end='\t')
				prev_night=mass
				print()
	print()

	plt.savefig('{}Overall_Mass.png'.format(storage_directory))
	#plt.show()

	print('Plotting Setpoints')
	plt.figure(figsize=(12,8))
	plt.title('Setpoints')
	for meta in setpoint_storage_dict:
		if meta !='time':
			print(meta)
			col=run_config_dict['Chunk Color'][meta]
			plt.plot(setpoint_storage_dict['time'],setpoint_storage_dict[meta],label=meta, color=col)
	plt.legend()
	for index_set in index_set_list:
		plt.axvspan(storage_dict['time'][index_set[0]],storage_dict['time'][index_set[1]], color='gray', alpha=1-light_list[index_set[0]]/run_config_dict['Light Regiment']['Day_Level'])
	#plt.gcf().tight_layout()
	ax=plt.gca()
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

	# Put a legend to the right of the current axis
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	plt.savefig('{}Setpoints.png'.format(storage_directory))

	print('Plotting Error')
	plt.figure(figsize=(12,8))
	plt.title('Error')
	for meta in chunk_track_dict:
		col=run_config_dict['Chunk Color'][meta]
		plt.plot(prod_error_storage_dict['time'],prod_error_storage_dict[meta],label=meta, color=col)
		plt.plot(deg_error_storage_dict['time'],deg_error_storage_dict[meta],'--', color=col)
	plt.legend()
	for index_set in index_set_list:
		plt.axvspan(storage_dict['time'][index_set[0]],storage_dict['time'][index_set[1]], color='gray', alpha=1-light_list[index_set[0]]/run_config_dict['Light Regiment']['Day_Level'])
	#plt.gcf().tight_layout()
	ax=plt.gca()
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

	# Put a legend to the right of the current axis
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	plt.savefig('{}Overall_Error.png'.format(storage_directory))

	print('Plotting Exchanges')
	plt.figure(figsize=(12,8))
	plt.title('Exchanges')
	for exchange in exchange_storage_dict:
		if exchange is not 'time':
		#col=run_config_dict['Chunk Color'][meta]
			plt.plot(exchange_storage_dict['time'], exchange_storage_dict[exchange], label=exchange)
		#plt.plot(exchange_storage_dict['time'],deg_error_storage_dict[meta],'--', color=col)
	plt.legend()
	for index_set in index_set_list:
		plt.axvspan(storage_dict['time'][index_set[0]],storage_dict['time'][index_set[1]], color='gray', alpha=1-light_list[index_set[0]]/run_config_dict['Light Regiment']['Day_Level'])
	#plt.gcf().tight_layout()
	ax=plt.gca()
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

	# Put a legend to the right of the current axis
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	plt.savefig('{}Overall_Exchanges.png'.format(storage_directory))
	#plt.show()

def reaction_transcript_investigator(fixed_chunked_model, reaction_id, reaction_gene_dict, max_flux, total_gene_set, conversion_file, gene_max_dict, storage_directory, model_parameter_file, run_config_dict):
	#this function interrogates a specific reaction to de-convolute the genes associated with the reaction
	print('Analysing {}'.format(reaction_id))
	bnds=fixed_chunked_model.reactions.get_by_id(reaction_id).bounds
	base_rule_string=reaction_gene_dict[reaction_id]['Rule']
	print()
	print('Gene Rule')
	print(base_rule_string)
	gene_time_value_dict={}
	lower_bound_list=[]
	upper_bound_list=[]

	critical_set=set([])

	x_space=np.linspace(0,24,240)
	for test_time in x_space:
		rule_string=copy.deepcopy(base_rule_string)
		gene_value_dict=gene_value_dict_maker(float(test_time), total_gene_set, conversion_file, gene_max_dict, model_parameter_file, run_config_dict)
		min_gene_value=1
		for gene in reaction_gene_dict[reaction_id]['Genes']:
			if gene not in gene_time_value_dict:
				gene_time_value_dict[gene]=[]
			gene_value=gene_value_dict.get(gene,1)
			min_gene_value=min(gene_value, min_gene_value)
			gene_time_value_dict[gene].append(gene_value)
			rule_string=rule_string.replace(gene,str(gene_value))
		#print()
		#print(rule_string)
		#print(eval(rule_string))
		#sys.exit()
		if rule_string!='':
			ret_val=min(max_flux*eval(rule_string),max_flux)
		else:
			ret_val=max_flux
		ret_val=max(ret_val,0)
		#new_bounds=[]
		if bnds[0]>0 or bnds[1]<0 or bnds==(0,0):
			new_bounds=bnds
		elif bnds[0]<0:
			new_bounds=(-ret_val,ret_val)
		else:
			new_bounds=(0,ret_val)
		#reaction.bounds=new_bounds
		lower_bound_list.append(new_bounds[0])
		upper_bound_list.append(new_bounds[1])

		#see if gene is critical at this time point:
		min_gene_value=max(0,min_gene_value)
		#print(min_gene_value)
		for test_gene in reaction_gene_dict[reaction_id]['Genes']:
			if gene_time_value_dict[test_gene][-1]<(min_gene_value+0.01):
				critical_set.add(test_gene)


	#for gene in reaction_gene_dict[reaction_id]['Genes']:
	#	if min(gene_time_value_dict[gene])<=0:
	#		critical_set.add(gene)
	#pprint(gene_time_value_dict.keys())
	pprint(critical_set)
	#print()
	#print(upper_bound_list)
	#sys.exit()
	work_fig=plt.figure()
	plt.axvspan(12,24, color='gray') #add darkness
	for transcript in critical_set:
		plt.figure(work_fig.number)
		plt.plot(x_space,gene_time_value_dict[transcript],label=transcript)
		plt.figure()
		plt.xlabel('Hours')
		plt.axvspan(12,24, color='gray', zorder=1) #add darkness
		plt.title(transcript)
		p=plt.plot(x_space,gene_time_value_dict[transcript],label='Fit')
		max_val=gene_max_dict[transcript]
		if transcript in conversion_file:
			transcript_list=conversion_file[transcript]
		else:
			transcript_list=[transcript]
		for transcript in transcript_list:
			gene='.'.join(transcript.split('.')[:2])
			#print(gene)
			#print(transcript)
			#print(transcript in gene_max_dict)
			#print(gene in model_parameter_file)
			plt.scatter(np.array(model_parameter_file[gene]['x'])%24, np.array(model_parameter_file[gene]['y'])/max_val, color=p[-1].get_color(), label='Data', zorder=5)

		plt.legend()
	plt.figure(work_fig.number)
	plt.title(reaction_id)
	plt.legend()
	ax2=plt.gca().twinx()
	col='xkcd:black'
	ax2.set_ylabel('Bounds')
	#ax2.plot(storage_dict['time'],total_mass_list,color=col,label='Total Mass')
	ax2.fill_between(x_space,upper_bound_list,lower_bound_list,alpha=0.2, label='Bounds')

	plt.savefig('{}{}.png'.format(storage_directory,reaction_id))


def total_chunk_mass_dump(run_config_dict, storage_dict, flux_export_dict, storage_directory, light_list):

	index_set_list=light_index(storage_dict,light_list)

	timepoint_index=0
	#print(len(storage_dict['time']))
	#print(len(flux_export_dict['Fluxes']))

	valid_index_list=[]
	valid_times=[]
	total_mass_dump_list=[]
	#sys.exit()
	#first check which reactions and indexes need to be examined
	while timepoint_index<len(storage_dict['time']):
		if isinstance(flux_export_dict['Fluxes'][timepoint_index], dict):
			valid_index_list.append(timepoint_index)
			valid_times.append(storage_dict['time'][timepoint_index])
			mass_dump=0
			flux_dict=flux_export_dict['Fluxes'][timepoint_index]
			for reaction in flux_dict:
				if '_chunk_export' in reaction:
					mass_dump+=flux_dict[reaction]
			total_mass_dump_list.append(mass_dump)

		timepoint_index+=1


	#pprint(limiting_reaction_set)
	#print(len(limiting_reaction_set))

	#plot all graphs
	print('Plotting Chunk Dumps')
	plt.figure()

	for index_set in index_set_list:
		plt.axvspan(storage_dict['time'][index_set[0]],storage_dict['time'][index_set[1]], color='gray', alpha=1-light_list[index_set[0]]/run_config_dict['Light Regiment']['Day_Level'], zorder=0)

	plt.plot(valid_times,total_mass_dump_list, zorder=10)
	#axis.set_ylabel(reaction_id,rotation=0, ha="left", rotation_mode="anchor")
	plt.ylabel('Total Mass Dumped')
	plt.xlabel('Time')
	plt.tight_layout()



def reaction_bound_plotter(run_config_dict, storage_dict, flux_export_dict, storage_directory, reaction_exclusion_list, light_list):

	index_set_list=light_index(storage_dict,light_list)

	timepoint_index=0
	#print(len(storage_dict['time']))
	#print(len(flux_export_dict['Fluxes']))
	limiting_reaction_set=set([])
	valid_index_list=[]
	valid_times=[]
	#sys.exit()
	#first check which reactions and indexes need to be examined
	while timepoint_index<len(storage_dict['time']):
		if isinstance(flux_export_dict['Fluxes'][timepoint_index], dict):
			valid_index_list.append(timepoint_index)
			valid_times.append(storage_dict['time'][timepoint_index])
			#print('*'*10)
			#print(storage_dict['time'][timepoint_index])
			flux_dict=flux_export_dict['Fluxes'][timepoint_index]
			bounds_dict=flux_export_dict['Bounds'][timepoint_index]
			for reaction_id in flux_export_dict['Fluxes'][timepoint_index]:
				if flux_dict[reaction_id]!=0 and flux_dict[reaction_id] in list(bounds_dict[reaction_id]):
					#print(flux_dict[reaction_id])
					#print(bounds_dict[reaction_id])
					#print('{}:({}, {}, {})'.format(reaction_id, bounds_dict[reaction_id][0], flux_dict[reaction_id], bounds_dict[reaction_id][1]))
					#print()
					#use both exact and wildcard logic for matching
					add_flag=True
					if reaction_id in reaction_exclusion_list:
						add_flag=False
					else:
						for entry in reaction_exclusion_list:
							if entry[0]=='*' and entry[1:] in reaction_id:
								add_flag=False
					if add_flag:
						limiting_reaction_set.add(reaction_id)
		timepoint_index+=1


	#pprint(limiting_reaction_set)
	#print(len(limiting_reaction_set))

	#plot all graphs
	print('Plotting restriction graphs')
	print("{} valid indexes".format(len(valid_index_list)))
	limiting_reaction_list=sorted(limiting_reaction_set)

	max_num=3
	list_of_list_of_reactions=[limiting_reaction_list[x:x+max_num] for x in range(0, len(limiting_reaction_list), max_num)]

	for sublist in list_of_list_of_reactions:
		fig, axes = plt.subplots(len(sublist), 1, sharex=True, squeeze=False)
		#pprint(axes)
		for index in range(len(sublist)):
			axis=axes[index][0]
			#plt.sca(axis)
			reaction_id=sublist[index]
			flux_list=[]
			lower_bound_list=[]
			upper_bound_list=[]
			for timepoint_index in valid_index_list:
				flux_list.append(flux_export_dict['Fluxes'][timepoint_index][reaction_id])
				lower_bound_list.append(flux_export_dict['Bounds'][timepoint_index][reaction_id][0])
				upper_bound_list.append(flux_export_dict['Bounds'][timepoint_index][reaction_id][1])

			flux_range=max(flux_list)-min(flux_list)
			for index_set in index_set_list:
				axis.axvspan(storage_dict['time'][index_set[0]],storage_dict['time'][index_set[1]], color='gray', alpha=1-light_list[index_set[0]]/run_config_dict['Light Regiment']['Day_Level'], zorder=0)

			axis.plot(valid_times,flux_list, zorder=10)
			axis.fill_between(valid_times,upper_bound_list,lower_bound_list,alpha=0.2, zorder=5)

			#axis.set_ylabel(reaction_id,rotation=0, ha="left", rotation_mode="anchor")
			axis.set_ylabel(reaction_id,rotation=0, va='top', ha='right')
			axis.yaxis.labelpad = 10
			axis.set_ylim(min(flux_list)-flux_range/2, max(flux_list)+flux_range/2)
			axis.set_xlim(min(valid_times),max(valid_times))
			#break

		fig.tight_layout()

	return limiting_reaction_list

def CNPlotter(storage_dict, chunk_output_dict):
	#print(storage_dict.keys())
	#print(storage_dict['time'])
	#print(np.zeros(len(storage_dict['time'])))
	C=np.zeros(len(storage_dict['time']))
	N=np.zeros(len(storage_dict['time']))
	for call_meta in storage_dict:
		if call_meta is not 'time':
			meta=call_meta+'_formation'
			C+=np.array(storage_dict[call_meta])*chunk_output_dict[meta]['Elements']['C']
			N+=np.array(storage_dict[call_meta])*chunk_output_dict[meta]['Elements']['N']
	#pprint(C)
	#pprint(N)
	plt.plot(storage_dict['time'],C/N)
	plt.title('C/N Ratio')
	plt.figure()
	#plt.show()

def setpoint_calculator(run_config_dict, sim_time, component_func_dict, DW_f):
	light_info=run_config_dict['Light Regiment']
	assert light_info['Light_Start']==0
	#assert light_info['Light_End']==12
	assert light_info['Day']==24
	#transform sim_time to measured time
	if sim_time%24<=12:
		transformed_time=sim_time%24
	else:
		transformed_time=-24+sim_time%24

	time_cell_DW=DW_f(transformed_time)
	setpoint_dict={}
	total_measured_mass=0
	unadjusted_remaining_mass=0
	for key in component_func_dict:
		if not isinstance(component_func_dict[key],float):
			#evaluate function at time point and add
			#print(type(component_func_dict[key]))
			setpoint_dict[key]=component_func_dict[key](transformed_time)/time_cell_DW
			total_measured_mass+=setpoint_dict[key]
		else:
			unadjusted_remaining_mass+=component_func_dict[key]

	#scale and add the remainders
	scaling_factor=(1-total_measured_mass)/unadjusted_remaining_mass
	for key in component_func_dict:
		if key not in setpoint_dict:
			setpoint_dict[key]=component_func_dict[key]*scaling_factor

	#pprint(component_func_dict)
	#pprint(setpoint_dict)
	#sys.exit()
	return setpoint_dict



def model_runner(run_config_dict, chunk_track_dict, chunk_error_integral_dict, chunked_model, fixed_chunked_model, total_gene_set, conversion_file, gene_max_dict, reaction_gene_dict, chunk_feasibilty_dict, ending_time, time_interval, storage_directory, max_flux, chunk_output_dict, model_parameter_file, component_func_dict, DW_f):

	#define the degradation mass-exchange minimization objective
	elm_mass_dict={'S':32.065, 'P':30.974, 'N':14.007, 'O':15.999, 'H':1.008, 'C':12.0111, 'Mg':24.305, 'Na':22.9897, 'Se':78.971, 'F':18.998, 'Fe':55.845}
	mass_dict=json.load(open('MetaboliteMassList.json'))
	mass_ex_objective=None
	for reaction in chunked_model.reactions:
		if len(reaction.metabolites)==1:
			#print(reaction)
			if 'CH_' in reaction.id or '_chunk_export' in reaction.id:
				factor=1
			else:
				factor=0
				assert len(reaction.metabolites)==1
				for metabolite in reaction.metabolites:
					if mass_dict[metabolite.id]!=0:
						factor=1000/mass_dict[metabolite.id]

			if mass_ex_objective==None:
				mass_ex_objective=factor*reaction.forward_variable
			else:
				mass_ex_objective+=factor*reaction.forward_variable
			mass_ex_objective+=factor*reaction.reverse_variable

	mass_ex_min_obj_obj=chunked_model.problem.Objective(mass_ex_objective, direction='min')


	flux_export_dict=None
	flux_export_flag=bool(run_config_dict['Export Flux Data'])
	if flux_export_flag:
		flux_export_dict={"Fluxes":[],"Bounds":[], "Times":[]}

	#sim start time
	#sim_time=23.9
	#sim_time=13
	sim_time=run_config_dict['Starting Time']
	setpoint_dict=setpoint_calculator(run_config_dict, sim_time, component_func_dict, DW_f)

	#initalize error and storage
	chunk_error_integral_dict, chunk_ratio_dict, chunk_production_error_dict, chunk_degradation_error_dict=error_calculator(chunk_track_dict, chunk_error_integral_dict, run_config_dict, setpoint_dict)

	storage_dict={'time':[]}
	for meta in chunk_track_dict:
		storage_dict[meta]=[]

	chunk_error_integral_storage_dict={'time':[]}
	for meta in chunk_error_integral_dict:
		chunk_error_integral_storage_dict[meta]=[]

	setpoint_storage_dict={'time':[]}
	for meta in setpoint_dict:
		setpoint_storage_dict[meta]=[]

	prod_error_storage_dict={'time':[]}
	for meta in chunk_track_dict:
		prod_error_storage_dict[meta]=[]

	deg_error_storage_dict={'time':[]}
	for meta in chunk_track_dict:
		deg_error_storage_dict[meta]=[]

	exchange_storage_dict={'time':[]}
	for exchange_reaction in run_config_dict['Tracked Exchange Reactions']:
		exchange_storage_dict[exchange_reaction]=[]

	light_list=[]
	total_mass_list=[]



	bound_reactions=set([])
	'''
	expression_list=[]
	expression_factors=[]
	for reaction in chunked_model.reactions:
		if len(reaction.metabolites)==1:
			if 'CH_' in reaction.id:
				print(reaction)
				print(reaction.flux_expression)
				print(type(reaction.flux_expression))
				expression_list.append(reaction.flux_expression)
				expression_factors.append(1)
			else:
				continue

			#print(reaction)
			#print(reaction.flux_expression)
			for meta in reaction.metabolites:
				form_dict=formula_dict_maker(meta.formula)
				if 'X' in form_dict or 'R' in form_dict:
					print(reaction)
					print(form_dict)
	print(expression_list)
	sys.exit()
	'''


	chunked_model, light_level=light_and_media_adjuster(run_config_dict, chunked_model, sim_time)

	print('Initalizing Run')
	ex_set=set([])
	while sim_time<ending_time:
		start_time=time.time()
		try:
		#if True:
			if not run_config_dict['FVA']:
				signal.alarm(100) #time out incase the finder hangs
			print('{:10.5f}'.format(sim_time),end='\t',flush=True)

			#track all values that need to be tracked
			storage_dict['time'].append(sim_time)
			light_list.append(light_level)
			total_mass=0
			for meta in chunk_track_dict:
				storage_dict[meta].append(chunk_ratio_dict[meta])
				total_mass+=chunk_track_dict[meta]
			total_mass_list.append(total_mass)

			chunk_error_integral_storage_dict['time'].append(sim_time)
			for meta in chunk_error_integral_dict:
				chunk_error_integral_storage_dict[meta].append(chunk_error_integral_dict[meta])

			setpoint_storage_dict['time'].append(sim_time)
			for meta in setpoint_dict:
				setpoint_storage_dict[meta].append(setpoint_dict[meta])

			prod_error_storage_dict['time'].append(sim_time)
			for meta in chunk_track_dict:
				prod_error_storage_dict[meta].append(chunk_production_error_dict[meta])

			deg_error_storage_dict['time'].append(sim_time)
			for meta in chunk_track_dict:
				deg_error_storage_dict[meta].append(chunk_degradation_error_dict[meta])

			chunked_model, light_level=light_and_media_adjuster(run_config_dict, chunked_model, sim_time)

			#chunked_model.reactions.CH_Energy_export.bounds=(0,1000)

			chunk_feasibilty_dict, chunk_track_dict, flux_export_dict, ex_rxn_dict=model_time_evaluator(run_config_dict, chunk_track_dict, chunk_production_error_dict, chunk_degradation_error_dict, chunked_model, fixed_chunked_model, total_gene_set, conversion_file, gene_max_dict, reaction_gene_dict, chunk_feasibilty_dict, sim_time, time_interval, flux_export_dict, mass_ex_min_obj_obj, mass_dict)


			#pprint(ex_rxn_dict)
			#print()
			for key in ex_rxn_dict:
				if abs(ex_rxn_dict[key])>1e-5:

					ex_set.add(key)
			#print()

			exchange_storage_dict['time'].append(sim_time)
			for exchange_reaction in run_config_dict['Tracked Exchange Reactions']:
				exchange_storage_dict[exchange_reaction].append(ex_rxn_dict.get(exchange_reaction,0))

			if 'EX_photonVis_e' not in ex_rxn_dict and ex_rxn_dict!={}:
				pprint(ex_rxn_dict)
				sys.exit()

			setpoint_dict=setpoint_calculator(run_config_dict, sim_time, component_func_dict, DW_f)

			chunk_error_integral_dict, chunk_ratio_dict, chunk_production_error_dict, chunk_degradation_error_dict=error_calculator(chunk_track_dict, chunk_error_integral_dict, run_config_dict, setpoint_dict)

			sim_time+=time_interval

			if round(sim_time,3)%1==0:
				print('savpoint',end='')
				json.dump(ex_rxn_dict,open('{}inprogress_flux_export_dict.json'.format(storage_directory),'w'), indent=4)

				json.dump(chunk_error_integral_dict,open('{}inprogress_error_integral.json'.format(storage_directory),'w'), indent=4)

				json.dump(setpoint_dict,open('{}inprogress_setpoint_storage.json'.format(storage_directory),'w'), indent=4)

				json.dump(chunk_ratio_dict,open('{}inprogress_ratio_dict.json'.format(storage_directory),'w'), indent=4)

				json.dump(total_mass_list,open('{}inprogress_total_mass_list.json'.format(storage_directory),'w'), indent=4)
			#print()
			#print('',end='\r')

		except CustomTimeExit:
			print()
			print('Failed, Leaving')
			sys.exit()
		end_time=time.time()
		print('e_t:{0:4.2f}'.format(end_time-start_time),end='\t')
		#print()
		print('',end='\n')
	signal.alarm(0)

	print('Exports')
	pprint(ex_set)
	#sys.exit()
	#pprint(flux_export_dict)

	print('Error Integrals')
	if False:
		for key in chunk_error_integral_storage_dict:
			print(key)
			print(chunk_error_integral_storage_dict[key][-3:])
			print()

	for key in chunk_error_integral_storage_dict:
		print('"{}":{},'.format(key,chunk_error_integral_storage_dict[key][-1]))

	print('\nStorage')
	if False:
		for key in storage_dict:
			print(key)
			print(storage_dict[key][-3:])
			print()

	for key in storage_dict:
		print('"{}":{},'.format(key,storage_dict[key][-1]))


	json.dump(flux_export_dict,open('{}flux_export_dict.json'.format(storage_directory),'w'), indent=4)

	json.dump(chunk_error_integral_storage_dict,open('{}error_integral.json'.format(storage_directory),'w'), indent=4)

	json.dump(setpoint_storage_dict,open('{}setpoint_storage.json'.format(storage_directory),'w'), indent=4)

	json.dump(storage_dict,open('{}storage_dict.json'.format(storage_directory),'w'), indent=4)

	json.dump(total_mass_list,open('{}total_mass_list.json'.format(storage_directory),'w'), indent=4)

	#sys.exit()
	if bool(run_config_dict['Plot C/N']):
		CNPlotter(storage_dict, chunk_output_dict)


	print()

	#plt.figure()


	if bool(run_config_dict['Plot Macro']):
		print('Plotting Mass')
		mass_and_error_plotter(run_config_dict, storage_dict, total_mass_list, light_list, prod_error_storage_dict, deg_error_storage_dict, storage_directory, exchange_storage_dict, setpoint_storage_dict)
		print('Plotting Chunk Dump')
		total_chunk_mass_dump(run_config_dict, storage_dict, flux_export_dict, storage_directory, light_list)
	if flux_export_flag and bool(run_config_dict['Plot Limits']):
		limiting_reaction_list=reaction_bound_plotter(run_config_dict, storage_dict, flux_export_dict, storage_directory, run_config_dict['Export Reaction Exclusion List'], light_list)

		limiting_reaction_list.extend(['PSIIblue','PSIIred','PSIblue','PSIred','RBPCh'])
		limiting_reaction_list=list(set(limiting_reaction_list))

		#plt.show()
		if bool(run_config_dict['Plot Transcripts']):
			for reaction_id in limiting_reaction_list:
				print(max_flux)
				if 'CH_' not in reaction_id or 'export' not in reaction_id:
					reaction_transcript_investigator(fixed_chunked_model, reaction_id, reaction_gene_dict, max_flux, total_gene_set, conversion_file, gene_max_dict, storage_directory, model_parameter_file, run_config_dict)

	if False:
		working_flux_dict=flux_export_dict['Fluxes'][55]
		key_list=sorted(working_flux_dict.items(), key=lambda x: x[1], reverse=True)
		for key in key_list:
			if abs(key[1])>1e-5:
				print(key)
				print(chunked_model.reactions.get_by_id(key[0]))
				print()
				print()


	if bool(run_config_dict['Plot Macro']) or bool(run_config_dict['Plot Limits']) or bool(run_config_dict['Plot C/N']) or run_config_dict['Plot Transcripts']:
		plt.show()
	sys.exit()


def formula_dict_maker(formula_string_in):
#built around 1 and 2 letter abbreviations
	#print(formula_string_in)
	formula_dict={}
	prev_is_alpha=None
	current_element=None
	current_amount=0
	for letter_info in enumerate(formula_string_in):
		#print(letter_info)
		cur_is_alpha=letter_info[1].isalpha()
		if cur_is_alpha:
			cur_is_caps=letter_info[1].isupper()
		else:
			cur_is_caps=None
		#print(prev_is_alpha)
		#print(cur_is_alpha)
		#print(cur_is_caps)
		if prev_is_alpha:
			if cur_is_alpha:
				if cur_is_caps:
					formula_dict[current_element]=float(1)
					current_element=letter_info[1]
				else:
					current_element='{}{}'.format(formula_string_in[letter_info[0]-1], letter_info[1])
					#print(current_element)
			else:
				current_amount+=str(letter_info[1])
		else:
			if cur_is_alpha:
				formula_dict[current_element]=float(current_amount)
				current_amount=''
				current_element=letter_info[1]
			else:
				current_amount+=str(letter_info[1])
		prev_is_alpha=cur_is_alpha
	if current_amount=='':
		current_amount=float(1)
	formula_dict[current_element]=float(current_amount)
	formula_dict.pop(None,None)
	#pprint(formula_dict)
	return formula_dict


if __name__ == "__main__":

	#extract models from the AICcModel module
	model_name_dict={
	'Kronecker_delta': 'Kronecker_delta',
 'constant_model': 'constant',
 'cos_mult_model': 'cos_mult_model',
 'cosine_1_term_model': 'cos1',
 'cosine_2_term_model': 'cos2',
 'fixed_decay': 'fixed_decay',
 'hat': 'hat',
 'hat_mult_model': 'hat_mult_model',
 'square_model': 'square_model',
 'square_width_model': 'square_width_model',
 'variable_decay': 'variable_decay'
 }

	#set plot font size
	plt.rcParams["font.family"] = "Times New Roman"
	matplotlib.rcParams.update({'font.size': 14})

	model_dict={}
	for i in dir(AICcModels):
		if i in model_name_dict:
			model_dict[model_name_dict[i]]=getattr(AICcModels, i)

	#pprint(model_dict)

	#change to change the config files
	run_config_dict=json.load(open('/home/homer/Documents/Transient_Modeling_Multiprocess_Moving_Setpoints/Wild_Type_Nh4_1_Ac_dir/Wild_Type_Nh4_1_Ac_Config.json','r'))

	storage_directory='{}/'.format(run_config_dict['Run Storage'])
	#make storage
	make_storage_directory(storage_directory)

	FW_data=pd.read_csv('Fresh_Weight.csv', index_col=0, dtype=np.float64) #data that's the mass of the cells at given time points
	#print(FW_data)
	DW_conversion_data=pd.read_csv('DW_Conversion_Data.csv', index_col=0, dtype=np.float64) #conversion data
	#print(DW_conversion_data)

	DW_data=pd.DataFrame(FW_data.loc[:,'fresh weight (pg/cell)']/DW_conversion_data.loc[:,'FW/DW Ratio'], index=FW_data.index, columns=['DW (pg/cell)'])
	#print(DW_data)
	DW_f=interp1d(DW_data.index,DW_data.loc[:,'DW (pg/cell)'])
	xspace=np.linspace(-12,12,1001)

	#process biomass data
	component_func_dict={}
	output_biomass_data_dict={}
	total_func_mass=np.zeros(len(xspace))
	for component in ['Protein','Starch','Chla','Chlb']:
		component_data=pd.read_csv('{}_Data.csv'.format(component), index_col=0, dtype=np.float64)
		component_func=interp1d(component_data.index, component_data.loc[:,'{} (pg/cell)'.format(component)])
		plt.plot(xspace,component_func(xspace)/DW_f(xspace),label=component)
		#output_biomass_data_dict['{} Mass (pg/cell)']=
		total_func_mass+=component_func(xspace)/DW_f(xspace)
		component_func_dict['CH_'+component]=component_func
	plt.plot(xspace,total_func_mass,label='Total')
	plt.legend()
	plt.savefig('Measured_Fraction_Plot.png')


	#sys.exit()
	#construct the initial values for the tracking dictionaries
	chunk_track_dict={}
	chunk_degradation_error_dict={}
	chunk_production_error_dict={}
	chunk_error_integral_dict={}
	total_goal_value=0
	for key in run_config_dict['Chunk Starting Information']:
		chunk_track_dict[key]=run_config_dict['Chunk Starting Information'][key]
		if 'Chunk Starting Integral' in run_config_dict and key in run_config_dict['Chunk Starting Integral']:
			chunk_error_integral_dict[key]=run_config_dict['Chunk Starting Integral'][key]
		else:
			chunk_error_integral_dict[key]=0
		if key not in component_func_dict:
			component_func_dict[key]=run_config_dict['Chunk Goal Information'][key]
		#total_goal_value+=run_config_dict['Chunk Goal Information'][key]

	#print('Goal Sum: {}'.format(total_goal_value))
	#mult_val=1/total_goal_value
	#for key in run_config_dict['Chunk Goal Information']:
	#	run_config_dict['Chunk Goal Information'][key]*=mult_val

	#new_goal_sum=0
	#for key in run_config_dict['Chunk Goal Information']:
	#	new_goal_sum+=run_config_dict['Chunk Goal Information'][key]

	#print('New Goal Sum: {}'.format(new_goal_sum))

	#sys.exit()

	#load or create chunked model
	model_name='iCre1355_auto.xml'
	chunked_model=model_chunker(model_name, run_config_dict)
	#chunked_model.solver='gurobi'
	#adjust model NGAM as needed
	chunked_model.reactions.get_by_id('ATPM_NGAM').bounds=(run_config_dict['NGAM'], run_config_dict['NGAM'])

	fixed_chunked_model=copy.deepcopy(chunked_model)
	#output model information
	chunk_output_dict={}
	#chunked_model=cobra.io.read_sbml_model('iCre1355_auto_chunked.xml')
	for reaction in chunked_model.reactions:
		if 'CH_' in reaction.id and ('_formation' in reaction.id or '_export' in reaction.id):
			chunk_output_dict[reaction.id]={'Name':reaction.name,'Metabolites':{}}
			for entry in reaction.metabolites:
				chunk_output_dict[reaction.id]['Metabolites'][entry.id]=reaction.metabolites[entry]
	json.dump(chunk_output_dict,open('chunk_used.json','w'),indent=4)
	chunk_output_dict=json.load(open('chunk_used.json'))
	for chunk in chunk_output_dict:
		if 'formation' in chunk:
			print(chunk)
			C=0
			N=0
			for metabolite_id in chunk_output_dict[chunk]['Metabolites']:
				if 'CH_' != metabolite_id[:3]:
					#print('{:25s} '.format(metabolite_id),end=' : ')
					#print(chunk_info[chunk]['Metabolites'][metabolite_id])
					form_dict=formula_dict_maker(chunked_model.metabolites.get_by_id(metabolite_id).formula)
					#pprint(form_dict)
					C-=chunk_output_dict[chunk]['Metabolites'][metabolite_id]*form_dict.get('C',0)
					N-=chunk_output_dict[chunk]['Metabolites'][metabolite_id]*form_dict.get('N',0)
					#print()
			print('C: {0:4.4f}\tN: {1:4.4f}'.format(C,N))
			chunk_output_dict[chunk]['Elements']={'C':C, 'N':N}
		#sys.exit()

	json.dump(chunk_output_dict,open('chunk_used.json','w'),indent=4)
	chunk_reaction_list=list(chunk_output_dict.keys())


	reaction_gene_dict={}
	#gene_value_dict={'Cre02.g120100.t1.2':10,'ChreCp049':5,'Cre02.g120150.t1.2':2}
	total_gene_set=set([])
	for reaction in chunked_model.reactions:
		gene_set=set([])
		for gene in reaction.genes:
			gene_set.add(gene.id)
			total_gene_set.add(gene.id)
		reaction_gene_dict[reaction.id]={'Genes':gene_set, 'Rule':gene_rule_dictionary_maker(reaction.gene_reaction_rule)}

	#code to update the values

	#print(run_config_dict['Transcript Source'])
	model_parameter_file=json.load(open(run_config_dict['Transcript Source'],'r'))
	conversion_file=json.load(open('conversion_table_organelle_gene_IDs.json','r')) #load for conversions between different gene ID schema

	print('Calculating Max For Genes')
	gene_max_dict={}
	x_list=np.linspace(0,24,200)
	for gene in total_gene_set:
		#print(gene)
		if gene in conversion_file:
			gene_key_list=conversion_file[gene]
			min_gene=np.inf
			for gene_key in gene_key_list:
				max_gene=0
				for x in x_list:
					gene_value=gene_value_finder(model_dict, model_parameter_file[gene_key],x)
					#print(gene_value)
					max_gene=max(max_gene,gene_value)
				min_gene=min(min_gene,max_gene)
			gene_max_dict[gene]=min_gene
		elif '.' in gene:
			gene_key='.'.join(gene.split('.')[:-2])
			#print(gene_key)
			#pprint(model_parameter_file[gene_key])
			max_gene=0
			for x in x_list:
				gene_value=gene_value_finder(model_dict, model_parameter_file[gene_key],x)
				#print(gene_value)
				max_gene=max(max_gene,gene_value)
			gene_max_dict[gene]=max_gene
	#pprint(gene_max_dict)
	#sys.exit()
	max_flux=run_config_dict['Max_flux']

	chunk_feasibilty_dict_name=run_config_dict['Chunk Feasibility File']
	chunk_feasibilty_dict={}
	try:
		chunk_feasibilty_dict=json.load( open('{}{}'.format(storage_directory, chunk_feasibilty_dict_name),'r'))
	except FileNotFoundError:
		print('Chunk feasibility not found at {}, creating new'.format(chunk_feasibilty_dict_name))
		chunk_feasibilty_dict={}




	print('Testing initial feasibility')
	base_chunk_feasibility_dict=json.load( open('{}{}'.format(storage_directory, run_config_dict['Base Chunk Feasibility File']),'r'))
	#chunk_feasibilty_dict=model_feasibility_tester(chunked_model,fixed_chunked_model, total_gene_set, conversion_file, gene_max_dict, reaction_gene_dict, chunk_feasibilty_dict, run_config_dict,time_set=[3.3000000000000003])

	#gets around a bug. Cobra hangs after a large number of evaluations. Splitting this into two helps for some reason.
	time_set_start=0
	interval=run_config_dict['Feasibility Save Interval']
	time_set_end=interval
	time_step=run_config_dict['Time Interval']
	while time_set_end<24:
		time_set_end=time_set_start+interval
		print('Running {} to {}'.format(time_set_start,time_set_end))
		time_set=np.arange(time_set_start,time_set_end,time_step)
		chunk_feasibilty_dict, overall_run_flag = model_feasibility_tester(chunked_model,fixed_chunked_model, total_gene_set, conversion_file, gene_max_dict, reaction_gene_dict, chunk_feasibilty_dict, run_config_dict, base_chunk_feasibility_dict, time_set)
		if overall_run_flag:
			json.dump(chunk_feasibilty_dict, open('{}{}'.format(storage_directory, chunk_feasibilty_dict_name),'w'), indent=4)
			chunk_feasibilty_dict=json.load( open('{}{}'.format(storage_directory, chunk_feasibilty_dict_name),'r'))
		time_set_start=time_set_end

	#time_set=np.linspace(12,24,121)
	#chunk_feasibilty_dict=model_feasibility_tester(chunked_model,fixed_chunked_model, total_gene_set, conversion_file, gene_max_dict, reaction_gene_dict, chunk_feasibilty_dict, run_config_dict, time_set)
	#json.dump(chunk_feasibilty_dict, open('{}{}'.format(storage_directory, chunk_feasibilty_dict_name),'w'), indent=4)



	if run_config_dict['FVA']:
		from cobra.flux_analysis import flux_variability_analysis
		obligatory_reactions={}
		#json.dump(obligatory_reactions,open('{}/FVA/obligatory_reactions.json'.format(storage_directory),'w'),indent=4)

	#pprint(chunk_error_integral_dict)

	model_runner(run_config_dict, chunk_track_dict, chunk_error_integral_dict, chunked_model, fixed_chunked_model, total_gene_set, conversion_file, gene_max_dict, reaction_gene_dict, chunk_feasibilty_dict, run_config_dict['Ending Time'],run_config_dict['Time Interval'], storage_directory, max_flux, chunk_output_dict, model_parameter_file, component_func_dict, DW_f)

	json.dump(run_config_dict,open('{}run_config_file.json'.format(storage_directory),'w'),indent=4)

	sys.exit()
