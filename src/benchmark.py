import numpy as np
import pandas as pd

from framed import FBA, solver_instance, Environment, essential_genes
from framed.cobra.ensemble import simulate_ensemble
from framed.solvers.solver import Status



def load_biolog_data(filename):
    return pd.read_csv(filename, sep='\t')


def simulate_biolog(model, medium, source, compounds, main_compounds, max_uptake=10, min_growth=0.1, verbose=False,
                    add_transporters=False, add_sinks=False, ensemble=False, voting_thresholds=None,
                    flavor=None):

    if ensemble:
        ensemble = model
        model = ensemble.model

    if flavor == 'seed':
        ex_rxn_format = 'EX_{}_e0'
    else:
        ex_rxn_format = 'R_EX_{}_e'

    env = Environment.from_compounds(medium, max_uptake=max_uptake, exchange_format='"{}"'.format(ex_rxn_format))
    constraints = env.apply(model, inplace=False, warning=False)

    for cmpd in main_compounds[source]:
        main_compound_rxn = ex_rxn_format.format(cmpd)
        constraints[main_compound_rxn] = (0, None)
    
    no_exchange = []
    not_in_model = []
    growth_pred = {}

    if flavor == 'seed':
        model_mets = {m_id[:-3] for m_id in model.metabolites}
    else:
        model_mets = {m_id[2:-2] for m_id in model.metabolites}

    new_rxns = []

    if add_transporters:
        model = model.copy()
        for met in compounds:
            met_e = ('{}_e0' if flavor == 'seed' else 'M_{}_e').format(met)
            met_c = ('{}_c0' if flavor == 'seed' else 'M_{}_c').format(met)
            if met_e not in model.metabolites and met_c in model.metabolites:
                if flavor == 'seed':
                    rxn_str = 'EX_{}: {} <-> [0, 0]'.format(met_e, met_c)
                else:
                    rxn_str = 'R_EX_{}_e: {} <-> [0, 0]'.format(met, met_c)
                new_rxns.append(model.add_reaction_from_str(rxn_str))

    if add_sinks:
        model = model.copy()
        for m_id in model.metabolites:
            if m_id.endswith('_c'):
                rxn_str = 'Sink_{}: {} --> '.format(m_id, m_id)
                model.add_reaction_from_str(rxn_str)

    solver = solver_instance(model)

    summary = {}

    for met in compounds:

        growth_pred[met] = [None]*len(voting_thresholds) if ensemble else None
        r_id = ex_rxn_format.format(met)

        if r_id in model.reactions:
            tmp = constraints[r_id] if r_id in constraints else (0, 0)
            constraints[r_id] = (-max_uptake, 0)
            if ensemble:
                growth_all = simulate_ensemble(ensemble, constraints=constraints, solver=solver, get_fluxes=False)
                growth_bool = [rate > min_growth if rate is not None else False for rate in growth_all]
                growth_pred[met] = [sum(growth_bool)/float(ensemble.size) >= t for t in voting_thresholds]
            else:
                sol = FBA(model, constraints=constraints, solver=solver)
                if sol.status == Status.OPTIMAL:
                    growth_pred[met] = sol.fobj > min_growth
                else:
                    growth_pred[met] = False
                constraints[r_id] = tmp

        else:
            if met in model_mets:
                no_exchange.append(met)
            else:
                not_in_model.append(met)

    if verbose and no_exchange:
        print 'No exchange reactions in model:', ' '.join(sorted(no_exchange))

    if verbose and not_in_model:
        print 'Metabolites not in model:', ' '.join(sorted(not_in_model))

    return growth_pred


def eval_result(entry, min_growth=0.1):
    if entry['predicted'] is not None and entry['predicted'] > min_growth:
        if entry['growth'] in {'++', '+'}:
            return 'TP'
        else:
            return 'FP'
    else:
        if entry['growth'] in {'--', '-'}:
            return 'TN'
        else:
            return 'FN'


def benchmark_biolog(model, medium, source, data, main_compounds, verbose=False, add_transporters=False, add_sinks=False,
                     ensemble=False, voting_thresholds=None, flavor=None, data_col='bigg_id'):

    compounds = data[data_col].dropna()
    growth = simulate_biolog(model, medium, source, compounds, verbose=verbose,
                             add_transporters=add_transporters, add_sinks=add_sinks,
                             ensemble=ensemble, voting_thresholds=voting_thresholds,
                             flavor=flavor, main_compounds=main_compounds)

    if ensemble:
        results = []
        for i in range(len(voting_thresholds)):
            growth_i = [(met, val[i]) for met, val in growth.items()]
            growth_df = pd.DataFrame(growth_i, columns=[data_col, 'predicted'])
            data_new = pd.merge(data, growth_df)
            data_new['result'] = data_new.apply(eval_result, axis=1)
            results.append(data_new)
        return results
    else:
        growth_df = pd.DataFrame(growth.items(), columns=[data_col, 'predicted'])
        data_new = pd.merge(data, growth_df)
        data_new['result'] = data_new.apply(eval_result, axis=1)
        return data_new


def summarize_benchmark(results, verbose=False):
    tp = results.count('TP')
    fn = results.count('FN')
    fp = results.count('FP')
    tn = results.count('TN')

    result = calculate_metrics(tp, tn, fp, fn)

    if verbose:
        print 'TP {:4d}  FN {:4d}'.format(tp, fn)
        print 'FP {:4d}  TN {:4d}'.format(fp, tn)
        print
        print 'Precision:   {:.4f}'.format(result['Precision'])
        print 'Accuracy:    {:.4f}'.format(result['Accuracy'])
        print 'Sensitivity: {:.4f}'.format(result['Sensitivity'])
        print 'Specificity: {:.4f}'.format(result['Specificity'])

    return result



def evaluate_reconstruction(real_reactions, new_reactions, universe, verbose=True):
    s_r = set(real_reactions)
    s_n = set(new_reactions)
    s_u = set(universe)

    tp = len(s_r & s_n)
    fn = len(s_r - s_n)
    fp = len(s_n - s_r)
    tn = len(s_u - (s_n | s_r))

    result = calculate_metrics(tp, tn, fp, fn)

    if verbose:
        print 'TP {:4d}  FN {:4d}'.format(tp, fn)
        print 'FP {:4d}  TN {:4d}'.format(fp, tn)
        print
        print 'Precision:   {:.4f}'.format(result['Precision'])
        print 'Accuracy:    {:.4f}'.format(result['Accuracy'])
        print 'Sensitivity: {:.4f}'.format(result['Sensitivity'])
        print 'Specificity: {:.4f}'.format(result['Specificity'])

    return result


def benchmark_essentiality(model, medium, in_vivo_essential, in_vivo_non_essential=None, verbose=False,
                           ensemble=False, voting_thresholds=None, flavor=None):

    if ensemble:
        ensemble = model
        model = ensemble.model

    if flavor == 'seed':
        ex_rxn_format = 'EX_{}_e0'
    else:
        ex_rxn_format = 'R_EX_{}_e'

    if medium is not None:
        env = Environment.from_compounds(medium, exchange_format='"{}"'.format(ex_rxn_format))
    else:
        env = Environment.complete(model)
        
    constraints = env.apply(model, inplace=False, warning=False)

    if ensemble:
        results = ensemble_essentiality(ensemble, constraints, voting_thresholds, min_growth=0.1)
        data = []
        for in_silico_essential in results:
            in_silico_non_essential = set(model.genes) - in_silico_essential
            res = essentialtity_eval(in_silico_essential, in_silico_non_essential, in_vivo_essential, in_vivo_non_essential, verbose)
            data.append(res)
        return data
    else:
        in_silico_essential = set(essential_genes(model, constraints=constraints, min_growth=0.1))
        in_silico_non_essential = set(model.genes) - in_silico_essential
        return essentialtity_eval(in_silico_essential, in_silico_non_essential, in_vivo_essential, in_vivo_non_essential, verbose)


def ensemble_essentiality(ensemble, constraints, voting_thresholds, min_growth=0.1):

        solver = solver_instance(ensemble.model)
        essential_all = []

        for i in range(ensemble.size):
            current = ensemble.get_constraints(i)

            if constraints:
                current.update(constraints)

            essential = essential_genes(ensemble.model, constraints=current, min_growth=min_growth, solver=solver)
            essential_all.extend(essential)


        results =[{gene for gene in ensemble.model.genes
                  if essential_all.count(gene)/float(ensemble.size) >= t}
                  for t in voting_thresholds]

        return results



def essentialtity_eval(in_silico_essential, in_silico_non_essential, in_vivo_essential, in_vivo_non_essential=None, verbose=False):

    in_vivo_essential = set(in_vivo_essential)

    if in_vivo_non_essential is None:
        model_genes = in_silico_essential | in_silico_non_essential
        in_vivo_non_essential = model_genes - in_vivo_essential

    tp = len(in_silico_essential & in_vivo_essential)
    fp = len(in_silico_essential & in_vivo_non_essential)
    fn = len(in_silico_non_essential & in_vivo_essential)
    tn = len(in_silico_non_essential & in_vivo_non_essential)

    result = calculate_metrics(tp, tn, fp, fn)

    result['TP'], result['TN'], result['FP'], result['FN'] = tp, tn, fp, fn

    if verbose:
        print 'TP {:4d}  FN {:4d}'.format(tp, fn)
        print 'FP {:4d}  TN {:4d}'.format(fp, tn)

    return result


def calculate_metrics(tp, tn, fp, fn):
    return {
    	'TP': tp, 'TN': tn, 'FP': fp, 'FN': fn,
        'Sensitivity': tp / float(tp + fn) if tp + fn > 0 else 0,
        'Specificity': tn / float(fp + tn) if fp + tn > 0 else 0,
        'Precision': tp / float(tp + fp) if tp + fp > 0 else 0,
        'Accuracy': (tp + tn) / float(tp + tn + fp + fn) if tp + tn + fp + fn > 0 else 0,
        'F1-score': 2 * tp / float(2 * tp + fp + fn) if tp + fp + fn > 0 else 0,
    }