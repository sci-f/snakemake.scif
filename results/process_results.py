#/usr/bin/env python3

import plotly.graph_objs as go
from glob import glob
import pandas
import json
import re
import os


df = pandas.DataFrame()
pwd = os.getcwd()
techs = ['singularity', 'shifter', 'charliecloud', 'docker', 'runc']

# Data goes into web directory
data_dir = '%s/web/data' %pwd
if not os.path.exists(data_dir):
    os.mkdir(data_dir)

columns = ['COMMAND',
           'ELAPSED_TIME_HMS',
           'AVERAGE_MEM',
           'FS_INPUTS',
           'MAX_RES_SIZE_KB',
           'FS_OUTPUTS',
           'PERC_CPU_ALLOCATED',
           'CPU_SECONDS_USED',
           'W_TIMES_SWAPPED',
           'SHARED_TEXT_KB',
           'ELAPSED_TIME_SECONDS',
           'NUMBER_SIGNALS_DELIVERED',
           'AVG_UNSHARED_STACK_SIZE',
           'SOCKET_MSG_RECEIVED',
           'SOCKET_MSG_SENT',
           'AVG_RESIDENT_SET_SIZE',
           'CONTEXT_SWITCHES']

for tech in techs:
    time_log = "%s/%s/snakemake-times.log" %(pwd, tech)
    rows = pandas.read_csv(time_log, sep='\t', header=None)
    rows.columns = columns
    rows.index = [tech]
    df = df.append(rows)


# A function to convert string H:M:S to minutes
def hms_to_minutes(time_string):
    # If there is a ., then we just have hours minutes
    time_string = time_string.split('.')[0]
    parts = time_string.split(':')
    if len(parts) == 2:
        m,s = parts
        h = 0
    else:
        h,m,s = parts
    h = int(h) * 60
    m = int(m)
    s = int(s) / 60.0
    return h+m+s

# A function to convert percent % to numerical
def perc_to_number(input_string):
    return [int(x.strip('%'))/100 for x in input_string]
    
df = df.rename(index=str, columns={"ELAPSED_TIME_SECONDS": "ELAPSED_TIME_SEC"})

# Make simple bar charts, with analyses sorted by their command
variables = df.columns[df.sum()!=0].tolist()

# Command isn't a variable :)
variables.pop(variables.index('COMMAND'))
variables.pop(variables.index('ELAPSED_TIME_HMS'))

# Show each container for each metric
for var in variables:
    traces = []
    subset = df[var]
    ydata = df[var].tolist()
    idx = df[var].index.tolist()
    if var == "ELAPSED_TIME_HMS":
        ydata = [hms_to_minutes(t) for t in ydata]
    if var == "PERC_CPU_ALLOCATED":
        ydata = perc_to_number(ydata)
    trace = go.Bar(x=idx,
                   y=ydata,
                   name=var)
    traces.append(trace)
    with open('%s/%s.json' %(data_dir,var),'w') as filey:
        json.dump(traces,filey)

# Finally, let's print an index.html file from the template!
template_file = "%s/web/template.html" %(pwd)
with open(template_file,'r') as filey:
    template = filey.read()


# These are info / about blocks to add to each, depending on the variable
varlookup = {'AVERAGE_MEM':'Maximum resident set size, or memory (RAM in KB)',
             'CONTEXT_SWITCHES':'Number of voluntary context-switches (e.g., while waiting for an I/O)',
             'CPU_SECONDS_USED':'Total number of CPU seconds used directly (in user mode, seconds)',
             'MAX_RES_SIZE_KB':"Maximum resident set size of the process (in Kilobytes)",
             'ELAPSED_TIME_SEC':'Elapsed real (wall clock) time used by the process (seconds)',
             'FS_OUTPUTS':'number of file system outputs',
             'FS_INPUTS':'number of file system inputs',
             'PERC_CPU_ALLOCATED':'Percentage of CPU allocated for'}

nodata = '<div class="heatmap-blank-metric"></div>'

# This is a header row of script names
data = '''<div class="heatmap-blank-metric"><a target="_blank" href="http://man7.org/linux/man-pages/man1/time.1.html">time</a></div>'''

variables.sort()

for v in range(len(variables)): # metric
    var = variables[v]
    data_path = '%s/web/data/%s.json' %(pwd, var)
    title = varlookup[var]
    newpath = '''<div class="heatmap-metric heatmap-level-%s" data-metric="%s" data-title="%s"><div class="heatmap-actual">%s</div></div>''' %(v+1,var,title,var)
    data = "%s\n%s" %(data,newpath)

template = template.replace('[[DATA]]',data)
index = "%s/web/index.html" %pwd
with open(index,'w') as filey:
    filey.writelines(template)


# Save final data file for (non empty) metrics
df.to_csv('%s/compiled_result.tsv' %pwd, sep='\t')
