import argparse
import os
import pandas as pd


##################################### arguments #####################################################
parser = argparse.ArgumentParser(description="script for write html report in TF ChIPseq analysis pipeline.\n",\
formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument("-sample_name", help = "sample_name")

parser.add_argument("-output_m_path", help = "output_m_path")

parser.add_argument("-GO_figure_suffix", help = "GO_figure_suffix, defalut=jpeg",  default="jpeg")

parser.add_argument("-ht", help = "high quality threshold using in selecting high quality peak, defalut=20",  default=20)

args = parser.parse_args()


########## parameters ##########
sample_name = args.sample_name
output_m_path = args.output_m_path
GO_figure_suffix = args.GO_figure_suffix
ht = args.ht # high quality threshold

### mapping result paramter
with open(os.path.join(output_m_path, "hisat2", "%s_alignmentstat"%sample_name),"r") as mapping_result:
    lines = mapping_result.readlines()
    raw_reads = lines[0].split(" ")[0]
    mapped_reads = lines[7].split(" ")[0]
    mapping_rate = lines[7].split(" (")[1].split(" ")[0]

### call peak result parameter
target_gene_file = pd.read_excel(os.path.join(output_m_path,"target_gene", "%s_target_gene.xlsx"%sample_name))
peak_summit_num = len(target_gene_file)
hq_peak_summit_num = len(target_gene_file[target_gene_file['signal']>ht])

### target gene file parameter
target_gene_file_sorted = pd.read_excel(os.path.join(output_m_path,"target_gene", "%s_target_gene_sorted.xlsx"%sample_name))
target_gene_num = len(set(target_gene_file_sorted['geneID']))
pr_target_gene_num = len(set(target_gene_file_sorted[target_gene_file_sorted['binding_type']=='promoter']['geneID'])) #promoter region target gene number
gb_targe_gene_num = len(set(target_gene_file_sorted[target_gene_file_sorted['binding_type']=='gene body']['geneID'])) # gene body target gene number

hq_target_gene_file_sorted = target_gene_file_sorted[target_gene_file_sorted['signal'] > ht].copy() # high quality

hq_target_gene_num = len(set(hq_target_gene_file_sorted['geneID']))
hq_pr_target_gene_num = len(set(hq_target_gene_file_sorted[hq_target_gene_file_sorted['binding_type']=='promoter']['geneID'])) # promoter region target gene number
hq_gb_targe_gene_num = len(set(hq_target_gene_file_sorted[hq_target_gene_file_sorted['binding_type']=='gene body']['geneID'])) # gene body target gene number


target_gene_table_html_head = """
<div id="target_gene_sorted">
<table>
    <thead>
        <tr>
            <th>chr</th>
            <th>summit</th>
            <th>binding_type</th>
            <th>signal</th>
            <th>geneID</th>
            <th>common_name</th>
            <th>description</th>
        </tr>
    </thead>
    <tbody>
"""
target_gene_table_html_body = ""
for i in range(len(target_gene_file_sorted)):
    chr = target_gene_file_sorted.iloc[i,0]
    summit = target_gene_file_sorted.iloc[i,1]
    binding_type = target_gene_file_sorted.iloc[i,2]
    signal = target_gene_file_sorted.iloc[i,3]
    geneID = target_gene_file_sorted.iloc[i,4]
    common_name = target_gene_file_sorted.iloc[i,5]
    description = target_gene_file_sorted.iloc[i,6]

    target_gene_table_html_body = target_gene_table_html_body + """<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>"""%(chr, summit, binding_type, signal, geneID, common_name, description)

target_gene_table_html_tail = """
    </tbody>
</table>
</div>
"""

target_gene_table_html = target_gene_table_html_head + target_gene_table_html_body + target_gene_table_html_tail

### GOresult parameter
figures = os.listdir(os.path.join(output_m_path, "GO"))
figures = [figure for figure in figures if GO_figure_suffix in figure]

# Initialize max values
max_C = -1
max_F = -1
max_P = -1

# Process each figure to extract category and number
for figure in figures:
    if "GO_result_C_" in figure:
        number = int(figure.split("C_")[1].split(".")[0])  # Extract number after C_
        max_C = max(max_C, number)
    elif "GO_result_F_" in figure:
        number = int(figure.split("F_")[1].split(".")[0])  # Extract number after F_
        max_F = max(max_F, number)
    elif "GO_result_P_" in figure:
        number = int(figure.split("P_")[1].split(".")[0])  # Extract number after P_
        max_P = max(max_P, number)


### write html content
content = """
<!DOCTYPE html>
<html>
  <head>
    <meta charset="utf-8">
    <title>%s</title>
    <link type="text/css" rel="stylesheet" href="../TF.css">
  </head>
  <body>
    <h1>%s</h1>
    <h2>data qality</h2>
    <p><a href="../fastqc/%s_fastqc.html">fastqc result</a></p>
    <h2>mapping reuslt</h2>
    <ul>
        <li>raw reads:%s</li>
        <li>mapped reads:%s</li>
        <li>mapping rate:%s</li>
    </ul>
    <a href="hisat2/%s_alignmentstat">mapping result download link </a>

    <h2> call peak result </h2>
    <ul>
        <li>The number of peak summit: %s</li>
        <li>The number of target gene:%s
            <ul>
                <li> binding in promoter region: %s</li>
                <li> binding in gene body: %s</li>
            </ul>
        </li>
        <br>
        <li>The number of high quality peak summit:%s</li>
        <li>The number of high quality target gene:%s</li>
            <ul>
                <li> binding in promoter region: %s</li>
                <li> binding in gene body: %s</li>
            </ul>
    </ul>

    <h2> target gene table </h2>
    %s

    <h2>GO Result</h2>

    <h2>BP</h2>
    <div id="GO_BP_result">
        <img src="GO/%s_GO_result_P_0.jpeg" class="GO_figure" id="bp_image">
    </div>
    <button id="bp_previous">Previous</button>
    <button id="bp_next">Next</button>
    <br>

    <h2>MF</h2>
    <div id="GO_MF_result">
        <img src="GO/%s_GO_result_F_0.jpeg" class="GO_figure" id="mf_image">
    </div>
    <button id="mf_previous">Previous</button>
    <button id="mf_next">Next</button>
    <br>

    <h2>CC</h2>
    <div id="GO_CC_result">
        <img src="GO/%s_GO_result_C_0.jpeg" class="GO_figure" id="cc_image">
    </div>
    <button id="cc_previous">Previous</button>
    <button id="cc_next">Next</button>
    <br>
    <p> all GO result in <a href="GO/%s_GO_result.xlsx">all go result download link </a> </p> 
    
    <script>
        // script to control GO figure
        const max_P = %s; // Example max value for P
        const max_F =%s; // Example max value for F
        const max_C = %s; // Example max value for C

        let current_P = 0;
        let current_F = 0;
        let current_C = 0;

        function updateButtons() {
            document.getElementById('bp_previous').disabled = current_P === 0;
            document.getElementById('bp_next').disabled = current_P === max_P;

            document.getElementById('mf_previous').disabled = current_F === 0;
            document.getElementById('mf_next').disabled = current_F === max_F;

            document.getElementById('cc_previous').disabled = current_C === 0;
            document.getElementById('cc_next').disabled = current_C === max_C;

            if (max_P === -1) document.getElementById('GO_BP_result').innerText = 'No result';
            if (max_F === -1) document.getElementById('GO_MF_result').innerText = 'No result';
            if (max_C === -1) document.getElementById('GO_CC_result').innerText = 'No result';
        }

        function showImages() {
            document.getElementById('bp_image').src = `GO/%s_GO_result_P_${current_P}.jpeg`;
            document.getElementById('mf_image').src = `GO/%s_GO_result_F_${current_F}.jpeg`;
            document.getElementById('cc_image').src = `GO/%s_GO_result_C_${current_C}.jpeg`;
            updateButtons();
        }

        document.getElementById('bp_next').onclick = function() {
            if (current_P < max_P) current_P++;
            showImages();
        };

        document.getElementById('bp_previous').onclick = function() {
            if (current_P > 0) current_P--;
            showImages();
        };

        document.getElementById('mf_next').onclick = function() {
            if (current_F < max_F) current_F++;
            showImages();
        };

        document.getElementById('mf_previous').onclick = function() {
            if (current_F > 0) current_F--;
            showImages();
        };

        document.getElementById('cc_next').onclick = function() {
            if (current_C < max_C) current_C++;
            showImages();
        };

        document.getElementById('cc_previous').onclick = function() {
            if (current_C > 0) current_C--;
            showImages();
        };

        // Initial display
        showImages();
    </script>

    <h2>motif analysis</h2>
    <p><a href="meme/meme_result/meme-chip.html">meme result</a></p>
</body>
</html>
"""%(sample_name, sample_name, sample_name, 
    raw_reads, mapped_reads, mapping_rate, sample_name, ##mapping result
    peak_summit_num, target_gene_num, pr_target_gene_num, gb_targe_gene_num, hq_peak_summit_num, hq_target_gene_num, hq_pr_target_gene_num, hq_gb_targe_gene_num, ## call peak result
    target_gene_table_html, ## target gene html
    sample_name, sample_name, sample_name, sample_name, max_P, max_F, max_C, sample_name, sample_name, sample_name,)


with open("%s/%s.html"%(output_m_path, sample_name),"w+") as f:
    f.write(content)