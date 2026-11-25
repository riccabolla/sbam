import math
import os
from datetime import datetime
from jinja2 import Environment, BaseLoader

HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>SBAM Assembly Report — {{ assembly_name }}</title>
    {% if not offline %}
    <script src="https://cdn.plot.ly/plotly-2.24.1.min.js"></script>
    {% endif %}
    <style>
        body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; 
        background: #f4f4f9; color: #333; margin: 0; padding: 20px; line-height: 1.5; }
        .container { max-width: 1260px; margin: 0 auto; background: white; padding: 30px; 
        border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }
        h1 { border-bottom: 2px solid #2c3e50; padding-bottom: 10px; color: #2c3e50; margin-bottom: 5px; }
        h2 { color: #34495e; margin-top: 40px; border-bottom: 1px solid #eee; padding-bottom: 5px; }
        .subtitle { color: #7f8c8d; font-size: 14px; margin-bottom: 30px; }
        
        .executive-summary { background: {{ result_color }}22; border: 2px solid {{ result_color }}; 
        padding: 20px; border-radius: 8px; margin-bottom: 30px; }
        .executive-summary h2 { margin: 0 0 10px 0; color: {{ result_color }}; border: none; padding: 0; }
        .executive-summary p { margin: 0; font-size: 15.5px; font-weight: 500; color: #2c3e50; }
        .notes-box { margin-top: 12px; padding: 10px 14px; background: rgba(255,255,255,0.75); border-left: 4px solid #f39c12; border-radius: 4px; font-size: 13.5px; color: #2c3e50; }
        
        .grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); gap: 20px; margin-bottom: 30px; }
        .card { padding: 20px; border-radius: 6px; background: #ecf0f1; border-left: 5px solid #3498db; }
        .card h3 { margin: 0 0 5px 0; font-size: 13px; text-transform: uppercase; color: #7f8c8d; }
        .card .value { font-size: 19px; font-weight: bold; color: #2c3e50; }
        .card .desc { font-size: 12px; color: #555; margin-top: 8px; }
        
        table { width: 100%; border-collapse: collapse; margin-top: 10px; font-size: 13px; }
        th, td { padding: 10px; text-align: left; border-bottom: 1px solid #ddd; }
        th { background-color: #f8f9fa; color: #2c3e50; }
        
        .badge-pass, .badge-acceptable, .badge-linear, .badge-rotated { background: #2ecc71; color: white; padding: 4px 8px; border-radius: 12px; 
        font-size: 11px; font-weight: bold; }
        .badge-atypical, .badge-unrotated, .badge-near_oric { background: #f39c12; color: white; padding: 4px 8px; border-radius: 12px; 
        font-size: 11px; font-weight: bold; }
        .badge-fail, .badge-warning { background: #e74c3c; color: white; padding: 4px 8px; border-radius: 12px; 
        font-size: 11px; font-weight: bold; }
        .badge-warn, .badge-no_data { background: #95a5a6; color: white; padding: 4px 8px; 
        border-radius: 12px; font-size: 11px; font-weight: bold; }
        
        .plot-container { background: #fff; border: 1px solid #ddd; border-radius: 8px; padding: 10px; margin-top: 20px; }

        .explanation-box-success { background: #eafaf1; border: 1px solid #2ecc71; border-radius: 6px; padding: 15px; margin-bottom: 20px; }
        .explanation-box-success h4 { margin: 0 0 10px 0; color: #27ae60; font-size: 14px; text-transform: uppercase; }
        .explanation-box-success p { margin: 0; font-size: 13px; color: #2c3e50; }

        .explanation-box-warning { background: #fef5e7; border: 1px solid #e67e22; border-radius: 6px; padding: 15px; margin-bottom: 20px; }
        .explanation-box-warning h4 { margin: 0 0 10px 0; color: #d35400; font-size: 14px; text-transform: uppercase; }
        .explanation-box-warning p { margin: 0; font-size: 13px; color: #2c3e50; }
        .footnote { font-size: 12px; color: #666; margin-top: 8px; font-style: italic; }
    </style>
    <script>
        function renderOfflinePolarSVG(containerId, thetaDeg, rawR, oricPos, terPos, seqLen) {
            var el = document.getElementById(containerId);
            if (!el) return;
            var size = 520, cx = size / 2, cy = size / 2, baseR = 155, amp = 55;
            var minR = Math.min.apply(null, rawR), maxR = Math.max.apply(null, rawR);
            var range = (maxR - minR) || 1;
            var pts = [];
            for (var i = 0; i < thetaDeg.length; i++) {
                var rad = (thetaDeg[i] - 90) * Math.PI / 180.0;
                var r = baseR + (((rawR[i] - minR) / range) - 0.5) * 2 * amp;
                var x = cx + r * Math.cos(rad), y = cy + r * Math.sin(rad);
                pts.push(x.toFixed(1) + "," + y.toFixed(1));
            }
            var oRad = ((oricPos / seqLen) * 360 - 90) * Math.PI / 180.0;
            var tRad = ((terPos / seqLen) * 360 - 90) * Math.PI / 180.0;
            var ox = cx + baseR * Math.cos(oRad), oy = cy + baseR * Math.sin(oRad);
            var tx = cx + baseR * Math.cos(tRad), ty = cy + baseR * Math.sin(tRad);

            el.innerHTML = '<svg viewBox="0 0 ' + size + ' ' + size + '" style="width:100%;max-width:520px;height:auto;display:block;margin:0 auto;">' +
                '<circle cx="' + cx + '" cy="' + cy + '" r="' + (baseR - amp) + '" fill="none" stroke="#eee" stroke-dasharray="3,3"/>' +
                '<circle cx="' + cx + '" cy="' + cy + '" r="' + (baseR + amp) + '" fill="none" stroke="#eee" stroke-dasharray="3,3"/>' +
                '<circle cx="' + cx + '" cy="' + cy + '" r="' + baseR + '" fill="none" stroke="#bdc3c7" stroke-width="5"/>' +
                '<polygon points="' + pts.join(" ") + '" fill="rgba(52,152,219,0.22)" stroke="#3498db" stroke-width="2"/>' +
                '<circle cx="' + ox.toFixed(1) + '" cy="' + oy.toFixed(1) + '" r="7" fill="#2ecc71" stroke="#fff" stroke-width="2"/>' +
                '<text x="' + (ox + 10).toFixed(1) + '" y="' + (oy - 8).toFixed(1) + '" font-size="13" font-weight="bold" fill="#27ae60">oriC (~' + Math.round(oricPos/1000) + ' kb)</text>' +
                '<circle cx="' + tx.toFixed(1) + '" cy="' + ty.toFixed(1) + '" r="7" fill="#e74c3c" stroke="#fff" stroke-width="2"/>' +
                '<text x="' + (tx + 10).toFixed(1) + '" y="' + (ty + 16).toFixed(1) + '" font-size="13" font-weight="bold" fill="#c0392b">ter (~' + Math.round(terPos/1000) + ' kb)</text>' +
                '<text x="' + cx + '" y="22" text-anchor="middle" font-size="12" fill="#7f8c8d">0 kb</text>' +
                '</svg>';
        }
    </script>
</head>
<body>
    <div class="container">
        <h1>SBAM Assembly Report</h1>
        <div class="subtitle">Assembly: <strong>{{ assembly_name }}</strong> | Generated: {{ date }}</div>
        
        <div class="executive-summary">
            <h2>Result: {{ result }}</h2>
            <p>{{ result_msg }}</p>
            {% if info_notes %}
            <div class="notes-box">
                <strong>Inspection Notes (Within Acceptance Thresholds):</strong>
                <ul style="margin: 5px 0 0 0; padding-left: 20px;">
                    {% for note in info_notes %}
                    <li>{{ note }}</li>
                    {% endfor %}
                </ul>
            </div>
            {% endif %}
            
            <h3 style="margin: 25px 0 10px 0; font-size: 15px; color: #2c3e50;">Assembly Overview ({{ num_contigs }} Contigs — All Screened for Circularity &amp; Base-Level Consensus)</h3>
            <table class="overview-table">
                <tr>
                    <th>Contig</th>
                    <th>Inferred Type</th>
                    <th>Length</th>
                    <th>Circularity</th>
                    <th>Base QC (&lt;50% Err | Masked)</th>
                    <th>Evaluation</th>
                </tr>
                {% for c in contig_summaries %}
                <tr>
                    <td><strong>{{ c.id }}</strong></td>
                    <td><span style="color: #2980b9; font-weight: 600;">{{ c.type }}</span></td>
                    <td>{{ "%.1f"|format(c.length_kb) }} kb</td>
                    <td><span class="badge-{{ c.j_stat|lower }}">{{ c.j_stat }}</span></td>
                    <td style="font-size: 12px;">{{ c.base_qc_str }}</td>
                    <td style="color: #555; font-weight: 500;">{{ c.evaluation }}</td>
                </tr>
                {% endfor %}
            </table>
        </div>
        
        <div class="grid">
            <div class="card" style="border-left-color: {% if j_status in ['PASS', 'LINEAR'] and failed_contigs_count == 0 %}#2ecc71{% elif j_status == 'NO_DATA' %}#95a5a6{% else %}#e74c3c{% endif %};">
                <h3>1. Seam &amp; Contig Continuity</h3>
                <div class="value">{{ primary_junc_score }} [{{ j_status }}]</div>
                <div class="desc">Measures read spanning across contig circularization boundaries ({{ failed_contigs_count }} failed contig(s) across assembly).</div>
            </div>
            <div class="card" style="border-left-color: {% if p_status == 'ACCEPTABLE' %}#2ecc71{% elif p_status == 'ATYPICAL' %}#f39c12{% elif p_status == 'WARNING' %}#e74c3c{% else %}#95a5a6{% endif %};">
                <h3>2. Replication Structure</h3>
                <div class="value">{{ primary_symmetry }} [{{ p_status }}]</div>
                <div class="desc">Measures smoothed bidirectional GC-skew symmetry and <i>dnaA</i> box enrichment over genomic background (&ge;500 kb).</div>
            </div>
            <div class="card" style="border-left-color: {% if not depth_sufficient_for_base_qc %}#95a5a6{% elif consensus_errors_per_mb <= 5.0 and low_concord_total_bp <= 1000 and masked_pct <= 5.0 %}#2ecc71{% else %}#f39c12{% endif %};">
                <h3>3. Primary Base Concordance</h3>
                <div class="value">QV &approx; {{ qv_display }} ({{ "%.2f"|format(consensus_errors_per_mb) }}/Mb)</div>
                <div class="desc">
                    Tested {{ "{:,}".format(evaluated_bases) }} bp on primary chromosome: <strong>{{ "{:,}".format(consensus_error_bases) }}</strong> sites &lt;50% support ({{ "%.2f"|format(consensus_errors_per_mb) }}/Mb) and <strong>{{ "{:,}".format(subconsensus_bias_bases) }}</strong> sub-consensus sites (50%–{{ "%.0f"|format(min_concordance * 100) }}%, {{ "%.1f"|format(subconsensus_bias_per_mb) }}/Mb).
                </div>
            </div>
            {% if has_genome_size %}
            <div class="card" style="border-left-color: {% if completeness_status == 'PASS' %}#2ecc71{% elif completeness_status == 'WARNING' %}#f39c12{% else %}#e74c3c{% endif %};">
                <h3>4. Size Completeness</h3>
                <div class="value">{{ "%.1f"|format(total_length_mb) }} Mb ({{ "%.1f"|format(completeness_pct) }}%)</div>
                <div class="desc">{{ completeness_msg }}</div>
            </div>
            {% endif %}
        </div>

        <h2>Circularity</h2>
        <table>
            <tr><th>Contig</th><th>Length (kb)</th><th>Depth (Copy Num)</th><th>Spanning / Clipped</th><th>Junction Score</th><th>Status</th></tr>
            {% for contig, metrics in junction.items() %}
            <tr>
                <td>{{ contig }}</td>
                <td>{{ "%.1f"|format(metrics.length / 1000) }}</td>
                <td>{{ metrics.avg_depth }}x <strong>({{ "%.1f"|format(metrics.copy_number) }}c)</strong></td>
                <td>{{ metrics.spanning_reads }} / {{ metrics.broken_reads }}</td>
                <td>{{ "%.2f"|format(metrics.spanning_score) }}</td>
                <td><span class="badge-{{ metrics.status|lower }}">{{ metrics.status }}</span></td>
            </tr>
            {% endfor %}
        </table>

        <h2>Replication Structure</h2>
        <table>
            <tr>
                <th>Contig</th>
                <th>Symmetry</th>
                <th>oriC Pos</th>
                <th>ter Pos</th>
                <th>DnaA Boxes (&plusmn;25 kb vs. Exp)</th>
                <th>oriC Offset (to 0 bp)</th>
                <th>Rotation</th>
                <th>Plausibility</th>
            </tr>
            {% for contig, metrics in physics.items() %}
            <tr>
                <td>{{ contig }}</td>
                <td>{{ metrics.symmetry }}°</td>
                <td>~{{ metrics.oric_pos }}</td>
                <td>~{{ metrics.ter_pos }}</td>
                <td>
                    {{ metrics.dnaa_boxes_50kb|default(0) }} / {{ "%.1f"|format(metrics.expected_dnaa_50kb|default(0.0)) }} exp
                    <strong>({{ "%.1f"|format(metrics.dnaa_fold_enrichment|default(0.0)) }}x, p={{ "%.3f"|format(metrics.dnaa_poisson_pval|default(1.0)) }})</strong>
                </td>
                <td>{{ "%.1f"|format(metrics.oric_offset_kb|default(0.0)) }} kb ({{ "%.1f"|format(metrics.oric_offset_deg|default(0.0)) }}°)</td>
                <td><span class="badge-{{ metrics.rotation_status|default('ROTATED')|lower }}">{{ metrics.rotation_status|default('ROTATED') }}</span></td>
                <td><span class="badge-{{ metrics.viability|lower }}">{{ metrics.viability }}</span></td>
            </tr>
            {% endfor %}
            {% if not physics %}
            <tr><td colspan="8" style="text-align: center;">No chromosomal contigs &ge;500 kb available for replication architecture analysis.</td></tr>
            {% endif %}
        </table>

        {% for contig, metrics in physics.items() %}
        {% if metrics.rotation_status in ['UNROTATED', 'NEAR_ORIC'] %}
        <div class="explanation-box-warning" style="margin-top: 15px;">
            <h4>Origin Rotation Notice (Contig {{ contig }}: {{ metrics.rotation_status }})</h4>
            <p>Predicted GC-skew minimum (<i>oriC</i> ~{{ metrics.oric_pos }} bp) is offset by <strong>{{ "%.1f"|format(metrics.oric_offset_kb) }} kb ({{ "%.1f"|format(metrics.oric_offset_deg) }}°)</strong> from coordinate 0 bp (informational; does not alter QC verdict). If submitting to databases that standardize chromosomes to start at <i>dnaA</i>, consider reorienting with <code>dnaapler</code> and re-running <code>sbam</code>.</p>
        </div>
        {% endif %}
        <div class="plot-container">
            <h3 style="text-align: center; color: #2c3e50; margin-bottom: 5px;">Genome Architecture Map: Contig {{ contig }}</h3>
            <div id="plot_{{ loop.index }}" style="width:100%; max-width:800px; height:560px; margin:0 auto;"></div>
            <script>
                (function() {
                    var plotId = "plot_{{ loop.index }}";
                    var raw_r = {{ metrics.plot_r | tojson }};
                    var theta = {{ metrics.plot_theta | tojson }};
                    var oric_pos = {{ metrics.oric_pos | tojson }};
                    var ter_pos = {{ metrics.ter_pos | tojson }};
                    var seq_len = {{ metrics.length | tojson }};

                    if (typeof Plotly === "undefined") {
                        renderOfflinePolarSVG(plotId, theta, raw_r, oric_pos, ter_pos, seq_len);
                        return;
                    }

                    var min_r = Math.min.apply(null, raw_r);
                    var max_r = Math.max.apply(null, raw_r);
                    var range = (max_r - min_r) || 1;
                    var baseline = 2.0;
                    var norm_r = raw_r.map(function(r) { return baseline + ((r - min_r) / range - 0.5) * 1.6; });

                    var contig_trace = { type: 'scatterpolar', r: Array(theta.length).fill(baseline), theta: theta, mode: 'lines', line: {color: '#bdc3c7', width: 6}, hoverinfo: 'none' };
                    var skew_trace = { type: 'scatterpolar', r: norm_r, theta: theta, mode: 'lines', line: {color: '#3498db', width: 2}, fill: 'tonext', fillcolor: 'rgba(52, 152, 219, 0.3)' };
                    var marker_trace = {
                        type: 'scatterpolar', r: [baseline, baseline], theta: [(oric_pos / seq_len) * 360, (ter_pos / seq_len) * 360], mode: 'markers+text',
                        marker: { color: ['#2ecc71', '#e74c3c'], size: 14, line: {color: 'white', width: 2} },
                        text: ['<b>oriC</b>', '<b>ter</b>'], textposition: ['top right', 'bottom left']
                    };
                    var layout = {
                        polar: { radialaxis: { visible: true, range: [0, 4] }, angularaxis: { direction: "clockwise" } },
                        showlegend: false, margin: { t: 40, b: 20, l: 40, r: 40 }
                    };
                    Plotly.newPlot(plotId, [contig_trace, skew_trace, marker_trace], layout, {responsive: true});
                })();
            </script>
        </div>
        {% endfor %}

        <h2 style="margin-top: 50px;">Base-Level Analysis (Primary Chromosome: {{ primary_contig_id }})</h2>
        <div class="{% if masked_pct > 5.0 or consensus_errors_per_mb > 5.0 or low_concord_total_bp > 1000 %}explanation-box-warning{% else %}explanation-box-success{% endif %}">
            <h4>Excluded Structural Regions ({{ "{:,}".format(masked_bases) }} bp / {{ "%.2f"|format(masked_pct) }}%) &amp; Discordance Breakdown</h4>
            <p style="margin-top: 8px;">
                Evaluated <strong>{{ "{:,}".format(evaluated_bases) }} bp ({{ "%.2f"|format(evaluated_pct) }}%)</strong> of primary contig <strong>{{ primary_contig_id }}</strong> against a global discordant baseline of <strong>{{ "{:,}".format(discordant_bases) }} sites ({{ "%.1f"|format(discordant_per_mb) }} / Mb, {{ "%.4f"|format(baseline_discordance_rate * 100) }}%)</strong>:
                <br>&bull; <strong>Majority Consensus Disagreement (&lt;50% read support):</strong> <strong>{{ "{:,}".format(consensus_error_bases) }}</strong> sites (<strong>{{ "%.2f"|format(consensus_errors_per_mb) }} / Mb</strong>; empirical concordance <strong>QV &approx; {{ qv_display }}</strong>) — candidate polishing errors, homopolymer indels, or uncollapsed repeat variants. <em>(Note: Read-concordance QV measures read-to-assembly agreement; errors shared by both reads and the assembly are invisible without an orthogonal reference.)</em>
                <br>&bull; <strong>Minority Sub-Consensus Read Disagreement (50%–{{ "%.0f"|format(min_concordance * 100) }}% read support):</strong> <strong>{{ "{:,}".format(subconsensus_bias_bases) }}</strong> sites (<strong>{{ "%.1f"|format(subconsensus_bias_per_mb) }} / Mb</strong>, {{ "%.4f"|format(subconsensus_bias_rate * 100) }}% of genome) — assembly base is supported by the read majority.
                {% if motif_count > 0 %}
                <br>&bull; <strong>Motif Family Coverage:</strong> {{ motif_count }} clustered motif family(ies) account for <strong>{{ "{:,}".format(explained_discordant_bases) }} / {{ "{:,}".format(discordant_bases) }} ({{ "%.1f"|format(explained_discordant_pct) }}%)</strong> of all discordant sites on the primary chromosome.
                {% endif %}
            </p>
        </div>
        
        {% if masked_regions %}
        <h3 style="font-size: 14px; color: #34495e; margin-top: 20px;">Largest Excluded Structural Regions (&ge; 50 bp)</h3>
        <table style="margin-bottom: 30px; font-size: 13px;">
            <tr><th>Region Start</th><th>Region End</th><th>Length (bp)</th><th>Mean Depth</th><th>Mean Concordance</th><th>Primary Classification</th></tr>
            {% for region in masked_regions %}
            <tr>
                <td>{{ region.start }}</td>
                <td>{{ region.end }}</td>
                <td><strong>{{ region.length }}</strong></td>
                <td>{{ "%.1f"|format(region.mean_depth|default(0.0)) }}x</td>
                <td>
                    {% if region.mean_concordance|default(100.0) < 90.0 %}
                    <strong style="color: #c0392b;">{{ "%.1f"|format(region.mean_concordance) }}%</strong>
                    {% else %}
                    {{ "%.1f"|format(region.mean_concordance|default(0.0)) }}%
                    {% endif %}
                </td>
                <td style="color: #e67e22;">{{ region.reason }}</td>
            </tr>
            {% endfor %}
        </table>
        {% endif %}
        
        <h3 style="font-size: 14px; color: #34495e; margin-top: 20px;">
            Clustered Motif Families (Primary Chromosome, Baseline {{ "%.4f"|format(baseline_discordance_rate * 100) }}%)
        </h3>
        <table>
            <tr>
                <th>Lead Canonical Motif</th>
                <th>Collapsed Family Members (Peak Member)</th>
                <th>Lead Hits / Total (Rate)</th>
                <th>Ref Orientation (+ / &minus;)</th>
                <th>Consensus Err (&lt;50%)</th>
                <th>Lead Fold (95% CI Lower)</th>
                <th>Assigned Sites (% Explained)*</th>
                <th>Orientation / Error Class</th>
                <th>Adj. p-value</th>
            </tr>
            {% for motif in motifs %}
            <tr>
                <td><strong>{{ motif.motif }}</strong> ({{ motif.kmer_size }}-mer)</td>
                <td style="font-size: 11.5px; color: #555;">
                    {% if motif.family_size|default(1) > 1 %}
                    {{ motif.family_size }} k-mers ({{ motif.related_motifs_str }})<br>
                    <em>Peak: {{ motif.peak_member_motif|default(motif.motif) }} ({{ "%.2f"|format(motif.max_member_error_rate|default(motif.error_rate) * 100) }}%, {{ "%.1f"|format(motif.max_member_fold|default(motif.fold_increase)) }}x)</em>
                    {% else %}
                    Single k-mer
                    {% endif %}
                </td>
                <td>{{ motif.discordant_count }} / {{ motif.occurrences }} ({{ "%.2f"|format(motif.error_rate * 100) }}%)</td>
                <td>
                    {% if motif.is_palindrome|default(false) %}
                    <em>Palindrome</em>
                    {% else %}
                    {{ motif.fwd_discordant|default(0) }} / {{ motif.rev_discordant|default(0) }}
                    {% endif %}
                </td>
                <td>
                    {% if motif.consensus_error_count|default(0) > 0 %}
                    <strong style="color: #c0392b;">{{ motif.consensus_error_count }}</strong>
                    {% else %}
                    <span style="color: #27ae60;">0</span>
                    {% endif %}
                </td>
                <td>{{ "%.1f"|format(motif.fold_increase) }}x (&ge;{{ "%.1f"|format(motif.ci_lower_fold|default(motif.fold_increase)) }}x)</td>
                <td>
                    <strong>{{ motif.assigned_unique_sites|default(motif.family_unique_sites) }}</strong> ({{ "%.1f"|format(motif.pct_discordant_explained|default(0.0)) }}%)
                    {% if motif.raw_family_sites is defined and motif.raw_family_sites != motif.assigned_unique_sites %}
                    <br><span style="font-size: 11px; color: #7f8c8d;">[{{ motif.raw_family_sites }} raw]</span>
                    {% endif %}
                </td>
                <td style="font-size: 11.5px;">
                    {% if motif.consensus_error_count|default(0) > 0 %}
                    <strong style="color: #c0392b;">{{ motif.artifact_class }}</strong>
                    {% elif motif.strand_biased|default(false) %}
                    <span style="color: #d35400; font-weight: 600;">{{ motif.artifact_class }}</span>
                    {% elif 'Undetermined' in motif.artifact_class|default('') %}
                    <span style="color: #7f8c8d;">{{ motif.artifact_class }}</span>
                    {% else %}
                    <span style="color: #2980b9;">{{ motif.artifact_class|default('Orientation-Balanced') }}</span>
                    {% endif %}
                </td>
                <td>{{ "%.2e"|format(motif.adj_p_value|default(motif.p_value)) }}</td>
            </tr>
            {% endfor %}
            {% if not motifs %}
            <tr><td colspan="9" style="text-align: center; color: #27ae60; font-weight: bold; padding: 20px;">
                {% if evaluated_bases > 0 %}
                No statistically enriched motif families (95% CI lower fold &ge; 2.0x) detected across {{ "{:,}".format(evaluated_bases) }} evaluated bases.
                {% else %}
                No unmasked bases were available for motif evaluation.
                {% endif %}
            </td></tr>
            {% endif %}
        </table>
        {% if motifs %}
        <div class="footnote">
            *<strong>Assigned Sites (% Explained):</strong> Discordant coordinates are partitioned mutually exclusively across families (assigned to the highest-ranked family) so the column sums to the exact deduplicated total ({{ "{:,}".format(explained_discordant_bases) }} / {{ "{:,}".format(discordant_bases) }} sites, {{ "%.1f"|format(explained_discordant_pct) }}%). Raw pre-deduplication family counts are shown in brackets when adjacent k-mer families overlap.<br>
            *<strong>Ref Orientation (+ / &minus;) &amp; Class:</strong> Counts discordant sites by whether the canonical motif appears on the forward (+) or reverse-complement (&minus;) reference strand. "Orientation-Asymmetric" requires two-sided binomial p &lt; 0.01 and &ge;80% single-orientation share; "Orientation-Balanced" requires n &ge; 20, p &ge; 0.05, and &le;65% share; smaller or intermediate counts are labeled "Undetermined".
        </div>
        {% endif %}
    </div>
</body>
</html>
"""


class DashboardBuilder:
    def __init__(self, outdir, assembly_path, offline=False):
        self.outdir = outdir
        self.assembly_name = os.path.basename(assembly_path)
        self.offline = offline
        self.last_verdict = "UNKNOWN"
        self.last_completeness_pct = None

    def generate_report(
        self,
        junction_metrics,
        physics_metrics,
        motif_results,
        masked_bases,
        masked_pct,
        masked_regions=None,
        genome_size=None,
        evaluated_bases=None,
        discordant_bases=None,
        consensus_error_bases=None,
        subconsensus_bias_bases=None,
        baseline_discordance_rate=None,
        consensus_error_rate=None,
        explained_discordant_bases=None,
        explained_discordant_pct=None,
        min_concordance=0.85,
        fidelity=None,
    ):
        print(" > [Report] Compiling HTML dashboard...")

        masked_regions = masked_regions or []
        motif_results = motif_results or []
        contig_base_metrics = {}

        if fidelity is not None:
            evaluated_bases = getattr(fidelity, "evaluated_bases", evaluated_bases)
            discordant_bases = getattr(fidelity, "discordant_bases", discordant_bases)
            consensus_error_bases = getattr(
                fidelity, "consensus_error_bases", consensus_error_bases
            )
            subconsensus_bias_bases = getattr(
                fidelity, "subconsensus_bias_bases", subconsensus_bias_bases
            )
            baseline_discordance_rate = getattr(
                fidelity, "baseline_discordance_rate", baseline_discordance_rate
            )
            consensus_error_rate = getattr(
                fidelity, "consensus_error_rate", consensus_error_rate
            )
            explained_discordant_bases = getattr(
                fidelity, "explained_discordant_bases", explained_discordant_bases
            )
            explained_discordant_pct = getattr(
                fidelity, "explained_discordant_pct", explained_discordant_pct
            )
            contig_base_metrics = getattr(fidelity, "contig_base_metrics", {}) or {}

        primary_contig_id = None
        if junction_metrics:
            primary_contig_id = max(
                junction_metrics.keys(), key=lambda k: junction_metrics[k]["length"]
            )

        primary_length = (
            junction_metrics[primary_contig_id]["length"]
            if primary_contig_id and primary_contig_id in junction_metrics
            else 0
        )
        primary_depth = (
            junction_metrics[primary_contig_id]["avg_depth"]
            if primary_contig_id and primary_contig_id in junction_metrics
            else 0.0
        )

        # Self-reconciling arithmetic
        if not evaluated_bases:
            evaluated_bases = max(0, primary_length - int(masked_bases or 0))
        if discordant_bases is None:
            discordant_bases = 0

        motif_consensus_sum = sum(
            m.get("consensus_error_count", 0) for m in motif_results
        )
        if consensus_error_bases is None or consensus_error_bases < motif_consensus_sum:
            consensus_error_bases = motif_consensus_sum

        subconsensus_bias_bases = max(
            0, int(discordant_bases) - int(consensus_error_bases)
        )

        if not baseline_discordance_rate and evaluated_bases > 0:
            baseline_discordance_rate = discordant_bases / evaluated_bases
        baseline_discordance_rate = float(baseline_discordance_rate or 0.0)

        if not consensus_error_rate and evaluated_bases > 0:
            consensus_error_rate = consensus_error_bases / evaluated_bases
        consensus_error_rate = float(consensus_error_rate or 0.0)

        subconsensus_bias_rate = (
            subconsensus_bias_bases / evaluated_bases if evaluated_bases > 0 else 0.0
        )

        eval_mb = evaluated_bases / 1_000_000 if evaluated_bases > 0 else 0.0
        consensus_errors_per_mb = (
            consensus_error_bases / eval_mb if eval_mb > 0 else 0.0
        )
        subconsensus_bias_per_mb = (
            subconsensus_bias_bases / eval_mb if eval_mb > 0 else 0.0
        )
        discordant_per_mb = discordant_bases / eval_mb if eval_mb > 0 else 0.0

        if evaluate_qv := evaluated_bases > 0:
            if consensus_error_bases <= 0:
                qv_val = round(-10.0 * math.log10(0.5 / evaluated_bases), 1)
                qv_display = f">{qv_val:.1f}"
            else:
                qv_val = round(
                    -10.0 * math.log10(consensus_error_bases / evaluated_bases), 1
                )
                qv_display = f"{qv_val:.1f}"
        else:
            qv_display = "N/A"

        if not explained_discordant_bases and motif_results:
            explained_discordant_bases = min(
                int(discordant_bases),
                sum(
                    m.get("assigned_unique_sites", m.get("family_unique_sites", 0))
                    for m in motif_results
                ),
            )
        explained_discordant_bases = int(explained_discordant_bases or 0)

        if not explained_discordant_pct and discordant_bases > 0:
            explained_discordant_pct = (
                explained_discordant_bases / discordant_bases
            ) * 100.0
        explained_discordant_pct = float(explained_discordant_pct or 0.0)

        low_concord_blocks = [
            r for r in masked_regions if r.get("mean_concordance", 100.0) < 90.0
        ]
        low_concord_block_count = len(low_concord_blocks)
        low_concord_total_bp = sum(r.get("length", 0) for r in low_concord_blocks)

        primary_junc_score = "N/A"
        j_status = "UNKNOWN"
        if primary_contig_id and primary_contig_id in junction_metrics:
            j_status = junction_metrics[primary_contig_id]["status"]
            primary_junc_score = (
                f"{junction_metrics[primary_contig_id]['spanning_score']:.2f}"
            )

        primary_symmetry = "N/A"
        p_status = "N/A"
        primary_rotation = "ROTATED"
        if primary_contig_id and primary_contig_id in physics_metrics:
            p_status = physics_metrics[primary_contig_id]["viability"]
            primary_symmetry = f"{physics_metrics[primary_contig_id]['symmetry']}°"
            primary_rotation = physics_metrics[primary_contig_id].get(
                "rotation_status", "ROTATED"
            )

        motif_count = len(motif_results)
        evaluated_pct = 100.0 - masked_pct
        depth_sufficient_for_base_qc = primary_depth >= 30.0

        failed_contigs_count = (
            sum(1 for m in junction_metrics.values() if m["status"] == "FAIL")
            if junction_metrics
            else 0
        )

        has_genome_size = genome_size is not None and genome_size > 0
        total_length = (
            sum(m["length"] for m in junction_metrics.values())
            if junction_metrics
            else 0
        )

        total_length_mb = total_length / 1_000_000
        genome_size_mb = (genome_size / 1_000_000) if has_genome_size else 0
        completeness_pct = (
            (total_length / genome_size * 100) if has_genome_size else 0.0
        )
        primary_pct = (
            (primary_length / genome_size * 100) if has_genome_size else 0.0
        )

        completeness_status = "UNKNOWN"
        completeness_msg = ""

        if has_genome_size:
            if completeness_pct < 90.0:
                completeness_status = "FAIL"
                completeness_msg = (
                    f"Assembly is missing >10% of expected sequence. "
                    f"Total size is {total_length_mb:.1f} Mb vs expected {genome_size_mb:.1f} Mb."
                )
            elif completeness_pct > 110.0:
                completeness_status = "WARNING"
                completeness_msg = (
                    "Total assembly is >10% larger than expected genome size "
                    "(may reflect accessory plasmids, uncollapsed repeats, or contamination)."
                )
            elif primary_pct < 85.0:
                completeness_status = "WARNING"
                completeness_msg = (
                    f"Total size is within range, but primary chromosome is only {primary_pct:.1f}% "
                    f"of expected size. Chromosome is likely fragmented."
                )
            else:
                completeness_status = "PASS"
                completeness_msg = (
                    "Total assembly size (90%–110%) and primary chromosome length (>=85%) match expectations."
                )

        num_contigs = len(junction_metrics) if junction_metrics else 0
        contig_summaries = []
        secondary_high_err_contigs = []

        if junction_metrics:
            for cid, j_metrics in junction_metrics.items():
                length_kb = j_metrics["length"] / 1000
                c_j_stat = j_metrics["status"]
                raw_class = j_metrics.get("classification", "Plasmid/Fragment")
                c_p_stat = physics_metrics.get(cid, {}).get("viability", "N/A")
                c_rot_stat = physics_metrics.get(cid, {}).get(
                    "rotation_status", "ROTATED"
                )

                c_base = contig_base_metrics.get(cid)
                if c_base:
                    c_err = c_base["consensus_error_bases"]
                    c_err_mb = c_base["consensus_errors_per_mb"]
                    c_mask_pct = c_base["masked_pct"]
                    base_qc_str = (
                        f"{c_err} err ({c_err_mb:.1f}/Mb) | {c_mask_pct:.2f}% masked"
                    )
                    if cid != primary_contig_id and c_err_mb > 5.0 and c_err >= 2:
                        secondary_high_err_contigs.append(
                            f"contig {cid} ({c_err} sites <50%, {c_err_mb:.1f}/Mb)"
                        )
                elif cid == primary_contig_id:
                    base_qc_str = f"{consensus_error_bases} err ({consensus_errors_per_mb:.2f}/Mb) | {masked_pct:.2f}% masked"
                else:
                    base_qc_str = "Not profiled"

                if c_j_stat in ("PASS", "LINEAR"):
                    if c_p_stat == "ACCEPTABLE":
                        evaluation = "Reliable Replicon Graph"
                        if c_rot_stat == "UNROTATED":
                            evaluation += " (Unrotated oriC >200 kb from 0 bp)"
                        elif c_rot_stat == "NEAR_ORIC":
                            evaluation += " (oriC within 200 kb of 0 bp)"
                    elif c_p_stat == "ATYPICAL":
                        evaluation = (
                            "Structurally Intact, Atypical replication symmetry"
                        )
                    elif c_p_stat == "WARNING":
                        evaluation = (
                            "Structurally Intact, Distorted replication symmetry"
                        )
                    elif "Sub-stoichiometric" in raw_class:
                        evaluation = "Structurally Intact (<0.7c: partial plasmid loss or extraction bias)"
                    else:
                        evaluation = "Structurally Intact"
                elif c_j_stat == "NO_DATA":
                    evaluation = "Insufficient Read Coverage (<5x)"
                else:
                    if "Misassembled" in raw_class:
                        evaluation = "Contradicted Junction"
                    elif "Debris" in raw_class:
                        evaluation = "Likely Low-Coverage Artifact"
                    else:
                        evaluation = "Unbridged / Linear Fragment"

                contig_summaries.append(
                    {
                        "id": cid,
                        "type": raw_class,
                        "length_kb": length_kb,
                        "j_stat": c_j_stat,
                        "base_qc_str": base_qc_str,
                        "evaluation": evaluation,
                    }
                )

        # Rate-Based Multi-Tier Verdict Logic + <30x Depth Guard
        info_notes = []

        if j_status == "NO_DATA":
            result = "INSUFFICIENT COVERAGE"
            result_color = "#7f8c8d"
            result_msg = (
                "Average read depth on the primary contig is below 5x. "
                "Cannot reliably evaluate junction spanning or base-level concordance."
            )
        elif j_status == "FAIL":
            result = "FAIL"
            result_color = "#c0392b"
            result_msg = (
                "Primary contig fails junction spanning thresholds. High probability of misassembly, "
                "incomplete circularization, or unresolved terminal repeats."
            )
        elif depth_sufficient_for_base_qc and consensus_errors_per_mb > 50.0:
            result = "FAIL (GROSS BASE DISCORDANCE)"
            result_color = "#c0392b"
            result_msg = (
                f"Primary chromosome exhibits a gross majority consensus error rate of "
                f"{consensus_errors_per_mb:.1f}/Mb ({consensus_error_bases} sites <50% support; QV ≈ {qv_display})."
            )
        else:
            if p_status == "ACCEPTABLE":
                result = "ACCEPT"
                result_color = "#27ae60"
                result_msg = (
                    f"Strong evidence of structural continuity and bidirectional replication architecture "
                    f"(empirical read-concordance QV ≈ {qv_display}, {consensus_errors_per_mb:.2f} consensus err/Mb)."
                )
            elif p_status in ("ATYPICAL", "WARNING"):
                result = "REVIEW"
                result_color = "#f39c12"
                result_msg = (
                    "Primary junction is intact, but replication symmetry deviates from standard "
                    "bidirectional geometry."
                )
            else:
                result = "ACCEPT WITH CAUTION"
                result_color = "#f39c12"
                result_msg = (
                    "Primary contig is structurally intact, but smaller than the 500 kb threshold "
                    "for chromosomal replication analysis."
                )

            review_triggers = []

            if not depth_sufficient_for_base_qc:
                if result == "ACCEPT":
                    result = "ACCEPT WITH CAUTION (LOW DEPTH <30x)"
                    result_color = "#f39c12"
                info_notes.append(
                    f"Primary depth ({primary_depth:.1f}x) is below 30x; insufficient depth for high-confidence base-level grading."
                )

            if failed_contigs_count > 0:
                review_triggers.append(
                    f"{failed_contigs_count} secondary contig(s) failed circularity"
                )
            if masked_pct > 5.0:
                review_triggers.append(
                    f"{masked_pct:.2f}% of the primary contig is structurally masked (>5% limit)"
                )

            if depth_sufficient_for_base_qc:
                if consensus_errors_per_mb > 5.0:
                    review_triggers.append(
                        f"consensus error rate is {consensus_errors_per_mb:.2f}/Mb "
                        f"({consensus_error_bases} sites <50% support, QV ≈ {qv_display}; threshold <=5.0/Mb)"
                    )
                elif consensus_error_bases > 0:
                    info_notes.append(
                        f"{consensus_error_bases} site(s) with <50% read support "
                        f"({consensus_errors_per_mb:.2f}/Mb, empirical QV ≈ {qv_display}) — within <=5.0/Mb acceptance threshold."
                    )

                if low_concord_total_bp > 1000 or low_concord_block_count >= 2:
                    worst_blk = min(
                        low_concord_blocks,
                        key=lambda b: b.get("mean_concordance", 100.0),
                    )
                    review_triggers.append(
                        f"{low_concord_block_count} excluded block(s) with <90% concordance totaling {low_concord_total_bp} bp "
                        f"(worst: {worst_blk['start']}–{worst_blk['end']} at {worst_blk['mean_concordance']}%)"
                    )
                elif low_concord_block_count == 1:
                    blk = low_concord_blocks[0]
                    info_notes.append(
                        f"1 contiguous excluded block ({blk['length']} bp at {blk['start']}–{blk['end']}) has "
                        f"{blk['mean_concordance']}% mean concordance ({blk['mean_depth']}x depth) — inspect in IGV for local repeat/variant."
                    )

                if secondary_high_err_contigs:
                    review_triggers.append(
                        "elevated consensus errors (>5/Mb) on plasmid/secondary "
                        + ", ".join(secondary_high_err_contigs)
                    )

            if motif_count > 0:
                info_notes.append(
                    f"{motif_count} motif family(ies) show sub-consensus read disagreement "
                    f"(explaining {explained_discordant_bases}/{discordant_bases} [{explained_discordant_pct:.1f}%] "
                    f"of sites <{min_concordance*100:.0f}% support)."
                )

            if primary_rotation == "UNROTATED":
                info_notes.append(
                    "Primary chromosome oriC is >200 kb from coordinate 0 bp (consider reorienting with dnaapler if required)."
                )

            if review_triggers:
                if result in ("ACCEPT", "ACCEPT WITH CAUTION", "ACCEPT (WITH NOTES)"):
                    result = "REVIEW"
                    result_color = "#f39c12"
                result_msg += " REVIEW FLAGGED: " + "; ".join(review_triggers) + "."
            elif result == "ACCEPT" and (
                consensus_error_bases > 0 or low_concord_block_count > 0
            ):
                result = "ACCEPT (WITH NOTES)"

        if has_genome_size and j_status != "NO_DATA":
            if completeness_status == "FAIL":
                result = "FAIL (INCOMPLETE)"
                result_color = "#c0392b"
                result_msg += (
                    f" Assembly is critically incomplete ({completeness_pct:.1f}% of expected size)."
                )
            elif completeness_status == "WARNING":
                if result in ("ACCEPT", "ACCEPT (WITH NOTES)", "ACCEPT WITH CAUTION"):
                    result = "REVIEW (SIZE ANOMALY)"
                    result_color = "#f39c12"
                if primary_pct < 85.0 and 90.0 <= completeness_pct <= 110.0:
                    result_msg += f" WARNING: Primary contig is only {primary_pct:.1f}% of expected genome size."
                else:
                    result_msg += f" WARNING: Total assembly size is {completeness_pct:.1f}% of expected genome size."

        self.last_verdict = result
        self.last_completeness_pct = (
            round(completeness_pct, 2) if has_genome_size else None
        )

        env = Environment(loader=BaseLoader(), autoescape=True)
        template = env.from_string(HTML_TEMPLATE)
        html_content = template.render(
            assembly_name=self.assembly_name,
            date=datetime.now().strftime("%Y-%m-%d %H:%M"),
            offline=self.offline,
            result=result,
            result_color=result_color,
            result_msg=result_msg,
            info_notes=info_notes,
            num_contigs=num_contigs,
            contig_summaries=contig_summaries,
            primary_contig_id=primary_contig_id,
            j_status=j_status,
            p_status=p_status,
            failed_contigs_count=failed_contigs_count,
            primary_junc_score=primary_junc_score,
            primary_symmetry=primary_symmetry,
            motif_count=motif_count,
            depth_sufficient_for_base_qc=depth_sufficient_for_base_qc,
            evaluated_bases=evaluated_bases,
            discordant_bases=discordant_bases,
            discordant_per_mb=discordant_per_mb,
            consensus_error_bases=consensus_error_bases,
            consensus_errors_per_mb=consensus_errors_per_mb,
            qv_display=qv_display,
            subconsensus_bias_bases=subconsensus_bias_bases,
            subconsensus_bias_per_mb=subconsensus_bias_per_mb,
            low_concord_block_count=low_concord_block_count,
            low_concord_total_bp=low_concord_total_bp,
            baseline_discordance_rate=baseline_discordance_rate,
            consensus_error_rate=consensus_error_rate,
            subconsensus_bias_rate=subconsensus_bias_rate,
            explained_discordant_bases=explained_discordant_bases,
            explained_discordant_pct=explained_discordant_pct,
            min_concordance=min_concordance,
            junction=junction_metrics,
            physics=physics_metrics,
            motifs=motif_results,
            masked_bases=masked_bases,
            masked_pct=masked_pct,
            evaluated_pct=evaluated_pct,
            masked_regions=masked_regions,
            has_genome_size=has_genome_size,
            total_length_mb=total_length_mb,
            completeness_pct=completeness_pct,
            completeness_status=completeness_status,
            completeness_msg=completeness_msg,
        )

        report_path = os.path.join(self.outdir, f"{self.assembly_name}_sbam_report.html")
        with open(report_path, "w", encoding="utf-8") as f:
            f.write(html_content)

        print(f" > [Report] Dashboard saved to: {report_path}")
        return report_path