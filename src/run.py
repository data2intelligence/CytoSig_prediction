#!/usr/bin/env python

import CytoSig
import pandas
import os, sys
import re
import numpy
import pathlib
import seaborn
import matplotlib.pyplot as plt

from scipy import stats
from glob import glob
from lifelines.fitters.coxph_fitter import CoxPHFitter
from lifelines import KaplanMeierFitter
from statsmodels.stats.multitest import multipletests
from scipy.cluster import hierarchy as hc
from matplotlib.backends.backend_pdf import PdfPages


###############################################################
# define paths
base_path = pathlib.Path(__file__).parent.absolute()
base_path = os.path.dirname(base_path)

data_path = os.path.join(base_path, 'data')
bulk_path = os.path.join(data_path, 'bulk')
sc_path = os.path.join(data_path, 'single_cell')
output_path = os.path.join(data_path, 'output')

###############################################################
# define parameters
alpha = 1E4
nrand = 1000
alternative='two-sided'
verbose_flag = False

figure_width = 7
font_size = 30

plt.rcParams.update({'font.size': font_size})

# default colors array
colors_default = plt.rcParams['axes.prop_cycle'].by_key()['color']

###############################################################
# load CytoSig signature
signature = os.path.join(sys.prefix, 'bin', 'signature.centroid')
signature = pandas.read_csv(signature, sep='\t', index_col=0)



###############################################################
# define functions

def category_background(expression, sep):
    """
    Create the background vector for single-cell analysis
    """
    
    if sep is None:
        background = expression.mean(axis=1)
    else:
        category = expression.groupby([v.split(sep)[0] for v in expression.columns], axis=1).apply(lambda v: v.mean(axis=1))
        background = category.mean(axis=1)
    
    background.name = 'background'
    return background



def heatmap(data, output, top_value=None, shrink=1, metric_type='correlation', flag_cluster=True):
    if type(data) == str: data = pandas.read_csv(data, sep='\t', index_col=0)
    
    if flag_cluster:
        data.fillna(0, inplace=True)
        z = hc.linkage(data, method='average', metric= metric_type)
        
        dendro = hc.dendrogram(z, labels=data.index.tolist(), no_plot=True)
        order_index = dendro['ivl']
        
        z = hc.linkage(data.transpose(), method='average', metric=metric_type)
        dendro = hc.dendrogram(z, labels= data.columns.tolist(), no_plot=True)
        order_columns = dendro['ivl']
        
        data = data.loc[order_index, order_columns]
    
    fig = plt.figure(figsize=data.shape, frameon=False)
    
    if top_value is None:
        g = seaborn.heatmap(data.transpose(), cmap="coolwarm", center=0, square=True, xticklabels=True, yticklabels=True, cbar_kws={"shrink": shrink})
    else:
        g = seaborn.heatmap(data.transpose(), cmap="coolwarm", center=0, vmax=top_value, vmin=-top_value, square=True, xticklabels=True, yticklabels=True, cbar_kws={"shrink": shrink})
                
    g.set_xticklabels(g.get_xmajorticklabels(), rotation=90)#, fontsize=font_size_label)
    g.set_yticklabels(g.get_ymajorticklabels(), rotation=0)#, fontsize=font_size_label)
    
    fig.savefig(output + '.pdf', bbox_inches='tight', transparent=True)
    #fig.savefig(output + '.svg', bbox_inches='tight', transparent=True)
    plt.close(fig)
    
    return data


def analyze_three_level_significance(out, levels, top_value = 4, cnt_thres = 5, qthres=0.05, vthres=1, merged_compare=True):
    """
    Statistical difference analysis between severe, mild, and healthy patient groups
    out: patient scores, such as CytoSig signalig activity or gene expression levels.
    top_value: the highest absolute z-score levels in plotting
    cnt_thres: minimum number of cells for each plotted sample
    qthres: threshold for significant corrected p-values
    vthres: threshold for minimum absolute z-scores to claim a significance event
    merged_compare: in the second comparison, whether merge severe and mild patients as one group against healthy controls
                    alternatively, just compare severe and healthy controls
    """
    
    result = pandas.read_csv(out, sep='\t', index_col=0)
    assert result.columns.value_counts().max() == 1

    # two types of comparison in analysis later
    label_comparison1 = '%s-%s' % (levels[0],  levels[1])
    
    if merged_compare:
        label_comparison2 = 'disease-%s' % levels[2]
    else:
        label_comparison2 = '%s-%s' % (levels[0],  levels[2])
    
    # (patient, cell type) group merge, merge on cell level
    flag = ['.'.join(v.split('.')[:2]) for v in result.columns]
    cntmap = pandas.Series(flag).value_counts()

    result = result.groupby(flag, axis=1).apply(lambda v: v.mean(axis=1))

    # only keep entities with sufficent number of cells
    result = result.loc[:, cntmap.index[cntmap >= cnt_thres]]

    result_group = result.groupby([v.split('.')[0] for v in result.columns], axis=1)

    writer = pandas.ExcelWriter(out + '.stat.xlsx', engine='xlsxwriter')
    pdf = PdfPages(out + '.heatmap.pdf')
    
    merge_signal_score = []
    
    mat_map_sig = {}
    cnt_map_sig = {}
    
    merge_severe_mild = []
    
    for cell_type, result in result_group:
        result.columns = [v.split('.')[1] for v in result.columns]
        result = result.loc[(result == 0).mean(axis=1) < 1]
        result = result.loc[result.std(axis=1) > 1e-3]
        
        # dataset specific fix
        if out.find('GSE145926') >= 0:
            result.columns = [v.replace('_C', '_') for v in result.columns]
        
        elif out.find('EGAS00001004571') >= 0:
            result.columns = [v.replace('_C19-CB-', '_') for v in result.columns]
        
        
        patient_group = result.groupby([v.split('_')[0] for v in result.columns], axis=1)

        if len( set(levels).intersection(patient_group.groups.keys())) != len(levels):
            # some patient group is missing for the current cell type
            print('jump %s by lacking of patient group' % cell_type)
            continue

        # order signal names 
        z = hc.linkage(result, method='average', metric='correlation')
        dendro = hc.dendrogram(z, labels=result.index.tolist(), no_plot=True)
        order_index = dendro['ivl']

        order_columns = []

        # order patients within each level group
        for level in levels:
            result_sub = patient_group.get_group(level)
            
            if result_sub.shape[1] > 1: # with something to order
                z = hc.linkage(result_sub.transpose(), method='average', metric='correlation')
                dendro = hc.dendrogram(z, labels= result_sub.columns.tolist(), no_plot=True)
                order_columns.extend(dendro['ivl'])
            else:
                order_columns.extend(result_sub.columns)

        result = result.loc[order_index, order_columns]

        fig = plt.figure(figsize=result.shape, frameon=False)
        g = seaborn.heatmap(result.transpose(), cmap="coolwarm", center=0, vmax=top_value, vmin=-top_value, square=True, xticklabels=True, yticklabels=True, cbar_kws={"shrink": 0.6})
        g.set_xticklabels(g.get_xmajorticklabels(), rotation=90)#, fontsize=font_size_label)
        g.set_yticklabels(g.get_ymajorticklabels(), rotation=0)#, fontsize=font_size_label)
        plt.title(cell_type, fontsize=font_size)

        pdf.savefig(fig, bbox_inches='tight', transparent=True)
        plt.close(fig)
             
        # differential analysis
        result_lst = []
        for level in levels: result_lst.append(patient_group.get_group(level))
        
        merge = []

        for gid in result_lst[0].index:
            arr_severe = result_lst[0].loc[gid]
            arr_mild = result_lst[1].loc[gid]
            arr_healthy = result_lst[2].loc[gid]
            
            if merged_compare:
                arr_disease = pandas.concat([arr_severe, arr_mild])
            else:
                arr_disease = arr_severe
            
            z_severe_mild, p_severe_mild = stats.ranksums(arr_severe, arr_mild)
            z_disease_healthy, p_disease_healthy = stats.ranksums(arr_disease, arr_healthy)
            
            merge.append(
                pandas.Series(
                    [arr_severe.median(), arr_mild.median(), arr_healthy.median(),
                     z_severe_mild, z_disease_healthy,
                     p_severe_mild, p_disease_healthy
                    ], name=gid,
                    index=['med.%s' % levels[0], 'med.%s' % levels[1], 'med.%s' % levels[2], 'z.%s' % label_comparison1, 'z.%s' % label_comparison2, 'p.%s' % label_comparison1, 'p.%s' % label_comparison2])
                )
        
        result_stat = pandas.concat(merge, axis=1, join='inner').transpose()
        
        signal_score = []
        
        for label_comparison in [label_comparison1, label_comparison2]:
            arr = result_stat['FDR.%s' % label_comparison] = multipletests(result_stat['p.%s' %label_comparison], method='fdr_bh')[1]
            signal_score.append(arr.mean())
            signal_score.append(arr.min())
        
        signal_score = pandas.Series(signal_score, index=[label_comparison1 + '.mean', label_comparison1 + '.min', label_comparison2 + '.mean', label_comparison2 + '.min'], name=cell_type)
        merge_signal_score.append(signal_score)
        
        result_stat.sort_values('p.%s' % label_comparison1, inplace=True)
        result_stat.to_excel(writer, sheet_name=cell_type[:31])
        
        # test for significance
        flag = (result_stat['FDR.%s' % label_comparison1] < qthres) & (result_stat['FDR.%s' % label_comparison2] < qthres)
        flag &= (
            (result_stat['med.%s' % levels[0]].abs() > vthres) | (result_stat['med.%s' % levels[1]].abs() > vthres) | (result_stat['med.%s' % levels[2]].abs() > vthres)
            )
        
        if sum(flag) > 0:
            print(cell_type)
            
            arr = result_stat.loc[:, 'z.%s' % label_comparison1]
            arr.name = cell_type
            merge_severe_mild.append(arr)
            
            result_stat = result_stat.loc[flag]
            result_stat.drop(['FDR.%s' % label_comparison1, 'FDR.%s' % label_comparison2], axis=1, inplace=True)
            result_stat.sort_values('z.%s' % label_comparison1, ascending=False).to_excel(writer, sheet_name= 'sig.' + cell_type)
            
            result_sub = result.loc[flag]
            result_sub.to_csv(out + '.' + cell_type, sep='\t', index_label=False)
            
            heatmap(result_sub, out + '.' + cell_type, top_value=top_value, shrink=0.5, flag_cluster=False)
            
            for gid in result_stat.index:
                prev = cnt_map_sig.get(gid, 0)
                cnt_map_sig[gid] = prev + 1
            
            mat_map_sig[cell_type] = result.groupby([v.split('_')[0] for v in result.columns], axis=1).median()
    
    pdf.close()
    

    # triplet pattern merge output
    if len(mat_map_sig) > 0:
        merge_severe_mild = pandas.concat(merge_severe_mild, axis=1, join='inner')
        
        if merge_severe_mild.shape[1] > 1:
            z = hc.linkage(merge_severe_mild.transpose(), method='average', metric='correlation')
            dendro = hc.dendrogram(z, labels= merge_severe_mild.columns.tolist(), no_plot=True)
            merge_severe_mild = merge_severe_mild.loc[:, dendro['ivl']]

        pandas.Series(cnt_map_sig).to_csv(out + '.cnt_map', sep='\t', index_label=False, header=False)
                
        merge = []
        
        for cell_type in merge_severe_mild.columns:
            triplet = mat_map_sig[cell_type]
            triplet.columns = cell_type  + '@' + triplet.columns
            merge.append(triplet.iloc[:, ::-1])
        
        pandas.concat(merge, axis=1, join='inner').to_csv(out + '.triplet', sep='\t', index_label=False)
    
    # output test statistics
    merge_signal_score = pandas.concat(merge_signal_score, axis=1, join='inner').transpose()
    merge_signal_score.sort_values(merge_signal_score.columns[0], inplace=True, ascending=False)
    merge_signal_score.to_csv(out + '.FDR', sep='\t', index_label=False)
    
    format_align = writer.book.add_format({'align': 'center'})
    format_number = writer.book.add_format({'num_format': '#,##0.000', 'align': 'center'})
    format_stat = writer.book.add_format({'num_format': '0.00E+00', 'align': 'center'})
        
    width = max(len('FDR.%s' % label_comparison1), len('FDR.%s' % label_comparison2))
        
    for worksheet in writer.sheets:
        worksheet = writer.sheets[worksheet]
            
        worksheet.set_column(0, result.shape[1], None, format_align)
        worksheet.set_column(1, 5, width, format_number)
        worksheet.set_column(6, 7, width, format_stat)
        worksheet.set_column(8, 9, width, format_number)
            
        worksheet.set_zoom(150)
    
    writer.close()



def survival_best_separation(data, pivot, z_continuous):
    """
    Get the signal cutoff for KM plot through the best-separation criteria.
    However, all statistical test will be based on continuous value without cutoffs..
    """
    cf = CoxPHFitter()
    
    arr = data[pivot]
    arr_rank = arr.sort_values(ascending=False)
    N = arr.shape[0]
    
    vthres_max = z_max = None
    
    margin = max(5, int(0.1*N))
    
    for i in range(margin-1, N-margin):
        vthres = arr_rank.iloc[i]
        
        data.loc[:, pivot] = (arr >= vthres)
        
        try:
            cf.fit(data, data.columns[0], event_col=data.columns[1])
        except:
            continue
        
        z = cf.summary.loc[pivot, 'z']
               
        if vthres_max is None or z * z_continuous > z_max * z_continuous:
            vthres_max = vthres
            z_max = z
    
    # recover the value        
    data.loc[:, pivot] = arr
    
    return vthres_max, z_max


def Survival_Analysis(clinical, data, signal, title, survival_type, output):
    cf = CoxPHFitter()
    
    data = pandas.read_csv(data, sep='\t', index_col=0)
    
    beta, se, zscore, pvalue = CytoSig.ridge_significance_test(signature, data, alpha, alternative, nrand, 1, True, False, verbose_flag)

    survival = pandas.read_csv(clinical, sep='\t', index_col=0)
    survival.index = survival.index.astype(str) # in case patient ID is integer
    
    mat = pandas.concat([survival, zscore.loc[signal]], axis=1, join='inner')

    cf.fit(mat, survival.columns[0], event_col=survival.columns[1])

    # wald test z-score
    z_Wald = cf.summary.loc[signal, 'z']

    # the CoxPH regression output two-sided p-values by default, we divide them by 2 to get one-sided
    p_oneside = cf.summary.loc[signal, 'p']/2

    # find the best separation cutoff only for KM plot, but not statistical test above
    thres, _ = survival_best_separation(mat, signal, z_Wald)

    flag = (mat[signal] > thres)

    fig = plt.figure(figsize=(figure_width, figure_width), frameon=False)
    kmf = KaplanMeierFitter()

    kmf.fit(mat.iloc[:,0].loc[flag], mat.iloc[:,1].loc[flag], label= 'High (n=%d)' % (sum(flag)))
    a1 = kmf.plot(ci_show=False, show_censors=True, color='red', linewidth=2)

    kmf.fit(mat.iloc[:,0].loc[~flag], mat.iloc[:,1].loc[~flag], label='Low (n=%d)' % (sum(~flag)))
    kmf.plot(ax=a1, ci_show=False, show_censors=True, color='blue', linestyle='--', linewidth=2)

    plt.xlabel(survival_type)
    plt.ylabel('Fraction')
    plt.legend(frameon=False)
    plt.title('p = %.2e' % p_oneside, fontsize=font_size)

    fig.savefig(output + '.pdf', bbox_inches='tight', transparent=True)
    plt.close(fig)
    

def Survival_Analysis_Tumor():
    """
    Analyze across tumor datasets
    """
    
    # data 1, Sunitinib multi-kinase inhibitor, targetting VEGF receptor, in metastatic clear-cell renal cell carcinoma (CCRCC)
    signal = 'VEGFA'
    title = 'Sunitinib CCRCC'
    survival_type = 'Progression-free survival (month)'

    fprefix = os.path.join(bulk_path, 'tumor', 'E-MTAB-3267')
    output = os.path.join(output_path, 'E-MTAB-3267.%s_%s' % (signal, title.replace(' ', '_')))

    clinical = fprefix + '.PFS'
    data = fprefix + '.norm_subtract.gz'

    Survival_Analysis(clinical, data, signal, title, survival_type, output)

    # data 2, Bevacizumab mono-therapy treatment in Glioblastoma (GBM)
    signal = 'VEGFA'
    title = 'Bevacizumab GBM'
    survival_type = 'Overall survival (month)'

    fprefix = os.path.join(bulk_path, 'tumor', 'GSE72951')
    output = os.path.join(output_path, 'GSE72951.%s_%s' % (signal, title.replace(' ', '_')))

    clinical = fprefix + '.OS.Bevacizumab'
    data = fprefix + '.self_subtract.gz'

    Survival_Analysis(clinical, data, signal, title, survival_type, output)

    # data 3, Atezolizumab anti-PDL1, requiring pre-existing interferon gamma signaling (IFNG), in metastatic urothelial carcinoma
    signal = 'IFNG'
    title = 'Atezolizumab Urothelial'
    survival_type = 'Overall survival (month)'

    fprefix = os.path.join(bulk_path, 'tumor', 'EGAS00001002556')
    output = os.path.join(output_path, 'EGAS00001002556.%s_%s' % (signal, title.replace(' ', '_')))

    clinical = fprefix + '.OS'
    data = fprefix + '.self_subtract.gz'

    Survival_Analysis(clinical, data, signal, title, survival_type, output)



def CytoSig_run_Inflam():
    """
    inflammatory diseases
    """

    # this is an output folder from the FDC framework, processed result files should look like *.diff.1
    fpath = os.path.join(bulk_path, 'inflam')
    
    # run CytoSig for all datasets
    diff_lst = glob(os.path.join(fpath, '*.diff.1'))
    
    for fprefix in diff_lst:
        print(os.path.basename(fprefix))
        
        data = pandas.read_csv(fprefix, sep='\t', index_col=0)
        beta, se, zscore, pvalue = CytoSig.ridge_significance_test(signature, data, alpha, alternative, nrand, 1, True, False, verbose_flag)
        zscore.to_csv(fprefix + '.signal', sep='\t', index_label=False)
        
        # if a dataset has individual sample information
        if os.path.exists(fprefix + '.sep.gz'):
            data = pandas.read_csv(fprefix + '.sep.gz', sep='\t', index_col=0)
            
            beta, se, zscore, pvalue = CytoSig.ridge_significance_test(signature, data, alpha, alternative, nrand, 1, True, False, verbose_flag)
            zscore.to_csv(fprefix + '.sep.signal', sep='\t', index_label=False)



def violinplot_one(handle, arr, i, col, alpha=0.6, flag_dot=False, alpha_dot=None, col_dot=None):
    parts = handle.violinplot(arr, [i], showmeans=False, showmedians=True, showextrema=False, widths=0.6)
    
    for pc in parts['bodies']:
        pc.set_facecolor(col)
        pc.set_edgecolor('black')
        pc.set_alpha(alpha)
    
    if 'cmeans' in parts:
        pc = parts['cmeans']
    elif 'cmedians' in parts:
        pc = parts['cmedians']
    else:
        pc = None
    
    if pc is not None:
        pc.set_linewidth(5)
        pc.set_edgecolor('black')
        pc.set_alpha(1)
        
    if flag_dot:
        if col_dot is None: col_dot = col
        if alpha_dot is None: alpha_dot = alpha
        
        x = numpy.random.normal(i, 0.1, size=arr.shape[0])
        handle.plot(x, arr, color=col_dot, marker='o', linestyle='none', markersize=5, alpha=alpha_dot)



def plot_canakinumab_GSE80060_response_activity():
    """
    systemic juvenile idiopathic arthritis (SJIA)
    """
    
    output = os.path.join(output_path, 'GSE80060.IL1B_Canakinumab_SJIA')
    
    conditions = ['100', '90', '70', '50', '30', '0', 'Placebo']
    
    # section 1: post versus pre activity
    result = os.path.join(bulk_path, 'inflam', 'GSE80060.MicroArray.GPL570.diff.1.sep.signal')
    result = pandas.read_csv(result, sep='\t', index_col=0)
    result = result.loc['IL1B']

    # standardize sample names
    info = pandas.DataFrame([v.split('@')[1].split('&')[0].split('_') for v in result.index], index=result.index)
    info.columns = ['Treatment', 'Response']
    
    flag_Canakinumab = (info.iloc[:,0] == 'Canakinumab') & (info.iloc[:,1] != 'nan')
    flag_Placebo = (info.iloc[:,0] == 'Placebo') & (info.iloc[:,1] == 'nan')
    info = info.loc[flag_Canakinumab | flag_Placebo]
    
    info.loc[flag_Canakinumab, 'Response'] = info.loc[flag_Canakinumab, 'Response'].apply(lambda v: v.replace('.0', ''))
    info = pandas.concat([info.loc[flag_Canakinumab, 'Response'], info.loc[flag_Placebo, 'Treatment']])
    
    result = result.loc[info.index]
    result_group = result.groupby(info)
    
    x = []
    y = []
    
    fig = plt.figure(figsize=(figure_width, figure_width), frameon=False)
    
    for i, conc in enumerate(conditions):
        result = result_group.get_group(conc)
        violinplot_one(plt, result, i, colors_default[0], flag_dot=True, alpha=0.5, alpha_dot=0.6)
        
        if conc == 'Placebo':
            conc = 0
        else:
            conc = int(conc)
        
        x.append(conc)
        y.append(result.median())
    
    r, p = stats.spearmanr(x, y)
    
    plt.title('r = %.2f, p = %.1e' % (r, p), fontsize=font_size)
    plt.axhline(0, linestyle='--', color='grey', lw=1)
    plt.xticks(range(len(conditions)), conditions, rotation=90)
    plt.tick_params(pad=10)
    
    plt.xlabel('Response (%)')
    plt.ylabel('IL1B activity diff day3')
    
    fig.savefig(output + '.pdf', bbox_inches='tight', transparent=True)
    plt.close(fig)

    
def plot_IFNK_GSE72754_response_activity():
    """
    systemic lupus erythematosus
    """
    
    output = os.path.join(output_path, 'GSE72754.IFN1_IFNK_SLE')
    
    # section 1: post versus pre activity
    result = os.path.join(bulk_path, 'inflam', 'GSE72754.MicroArray.HG-U133_Plus_2.diff.1.sep.signal')
    result = pandas.read_csv(result, sep='\t', index_col=0)
    result = result.loc['IFN1']
    
    info = pandas.DataFrame([v.split('&')[1].split()[0].split('_') for v in result.index], index=result.index, columns=['titer', 'month']).astype(int)
    
    r, p = stats.spearmanr(info['titer'], result)

    fig = plt.figure(figsize=(figure_width, figure_width), frameon=False)
    
    plt.plot(info['titer'], result, 'o', markersize=10)
    plt.title('r = %.2f, p = %.1e' % (r, p), fontsize=font_size)
    plt.xlabel('anti-IFNA Titer')
    plt.ylabel('IFN1 activity diff')
    
    plt.axhline(0, linestyle='--', color='grey')
    plt.axvline(100, linestyle='--', color='grey')
    
    plt.gca().ticklabel_format(axis="x", style="sci", scilimits=(0,0))
    
    fig.savefig(output + '.pdf', bbox_inches='tight', transparent=True)
    plt.close(fig)
    


def CytoSig_run_COVID19_singlecell():
    run_lst = [
        ['GSE145926', ['severe', 'mild', 'healthy']],
        ['EGAS00001004571_cohort1', ['severe', 'mild', 'control']],
    ]
    
    for dataset, levels in run_lst:
        data = os.path.join(sc_path, dataset + '.pickle.gz')
        output = os.path.join(output_path, dataset + '.signal')
        
        data = pandas.read_pickle(data)
        
        # normalize by health controls
        background = data.loc[:, [v.find(levels[2]) >= 0 for v in data.columns]]
        background = category_background(background, '.')
        data = data.subtract(background, axis=0)
        
        beta, se, zscore, pvalue = CytoSig.ridge_significance_test(signature, data, alpha, alternative, nrand, 1, True, False, verbose_flag)
        zscore.to_csv(output, sep='\t', index_label=False)
        
        qthres = 0.05
        # please see Methods section of our paper for neutrophils
        if dataset == 'EGAS00001004571_cohort1': qthres = 0.052
        
        analyze_three_level_significance(output, levels, qthres=qthres)




def analyze_COVID19_SC_Severity_joint_heatmap(top_value = 4):
    datasets = ['GSE145926', 'EGAS00001004571_cohort1']
    
    merge = []

    included = set()
    
    for dataset in datasets:
        output = os.path.join(output_path, dataset + '.signal')
        
        mat = pandas.read_csv(output + '.triplet', sep='\t', index_col=0)
        merge.append(mat)
        
        N = len(set([v.split('@',1)[0] for v in mat.columns]))
        
        cnt_map = pandas.read_csv(output + '.cnt_map', sep='\t', index_col=0, header=None)
        cnt_map = cnt_map.iloc[:, 0]
        cnt_map = cnt_map.loc[cnt_map >= N/3.0]
        
        included.update(cnt_map.index)
    
    merge = pandas.concat(merge, axis=1, join='inner')
    merge = merge.loc[included].transpose()
    
    post_s = re.compile('s$')
    
    lst = []
    for i, v in enumerate(merge.index):
        if (i-1) % 3 == 0:
            v = re.sub(post_s, '', v.split('@')[0])
        else:
            v = ''
        
        lst.append(v)
    
    merge.index = lst
    
    g = seaborn.clustermap(merge, cmap="coolwarm", row_cluster=False, metric='correlation', vmax=top_value, vmin=-top_value, center=0, xticklabels=True, yticklabels=True, cbar_kws={'label': 'Median Activity', "shrink": .5},
                               figsize=(1.1*merge.shape[1], 0.5*merge.shape[0]))#, dendrogram_ratio=(0.5, 0.1))
    
    plt.tick_params(pad=10)
    
    xlabel_lst = g.ax_heatmap.get_xticklabels()
    
    plt.setp(xlabel_lst, rotation=90)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
    
    merge = merge.loc[:, [v.get_text() for v in xlabel_lst]]
            
    for a in g.ax_row_dendrogram.collections: a.set_linewidth(2)
    for a in g.ax_col_dendrogram.collections: a.set_linewidth(2)
    
    out = os.path.join(output_path, 'COVID19.signal.triplet')
    
    merge.to_csv(out + '.sourcedata.csv')
    
    plt.savefig(out + '.significant.pdf', bbox_inches='tight', transparent=True)    


  

def main():
    # Figure 4e & Extended Data Figure 4 
    Survival_Analysis_Tumor()
    
    # prepare for Figure 4a-c
    CytoSig_run_Inflam()
    
    # Figure 4a
    plot_canakinumab_GSE80060_response_activity()
    
    # Figure 4b
    plot_IFNK_GSE72754_response_activity()
    
    # prepare for Figure 6 and generate Figure 6f
    CytoSig_run_COVID19_singlecell()
    
    # Figure 6g
    analyze_COVID19_SC_Severity_joint_heatmap()
    
    return 0

if __name__ == '__main__': main()
