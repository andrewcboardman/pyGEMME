import matplotlib.pyplot as plt
import seaborn as sns
from Bio.Phylo import draw

def plot_trace(trace):
    """Line plot of evolutionary conservation"""
    fig, ax = plt.subplots()
    plt.plot(trace.trace.values)
    plt.xlabel('Residue')
    plt.ylabel('Evolutionary trace')
    return ax

def plot_d_evol(d_evol):
    """Histogram of evolutionary distances to query"""
    fig, ax = plt.subplots()
    plt.hist(d_evol.d_evol.values, bins=100)
    plt.xlabel('Evolutionary distance')
    plt.ylabel('Number of sequences')
    return ax

def plot_fitness(fitness):
    """Heatmap of fitness scores"""
    fig, ax = plt.subplots()
    fitness = fitness.pivot(index='pos', columns='aa_mut', values='combi')
    sns.heatmap(fitness, cmap='viridis_r')
    plt.xlabel('Residue index')
    plt.ylabel('Substitution')
    #plt.colorbar()
    return ax

def make_plots(trace, d_evol, fitness, output_dir):
    p1 = plot_trace(trace)
    plt.savefig(f'{output_dir}/trace.svg')
    
    p2 = plot_d_evol(d_evol)
    plt.savefig(f'{output_dir}/d_evol.svg')

    p3 = plot_fitness(fitness)
    plt.savefig(f'{output_dir}/fitness.svg')
