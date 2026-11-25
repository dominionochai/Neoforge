/* NEOFORGE research explorer — vanilla JS, no build step, no framework.
   Data contract: frontend/data/case_<id>.json and benchmark_tesla.json,
   produced by `python -m neoforge.export_ui`. Those files are projections of
   content-addressed artifacts (every block carries the upstream artifact id +
   SHA-256 in its `source` field); this app renders them and performs no
   scientific computation of its own. The single interactive recomputation —
   the priority-rule playground — restates the published artifact rule and
   self-checks against the artifact's tier counts on every render. */
"use strict";

/* ---------------- tiny helpers ---------------- */

const $ = (sel, el) => (el || document).querySelector(sel);
const app = $("#app");

const ESC_MAP = { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" };
const esc = (s) => String(s ?? "").replace(/[&<>"']/g, (c) => ESC_MAP[c]);
const nfmt = (n) => new Intl.NumberFormat("en-GB").format(n);
const num = (x, d = 2) =>
  x === null || x === undefined ? "—" :
  (typeof x === "number" && !Number.isInteger(x) ? x.toFixed(d) : nfmt(x));
const shortSha = (s) => s ? s.slice(0, 12) + "…" : "—";

/** column-major artifact table -> array of row objects */
function rowsOf(tbl) {
  if (!tbl || !tbl.columns || !tbl.columns.length) return [];
  const cols = tbl.columns.map((c) => tbl[c]);
  const n = cols[0].length;
  const out = new Array(n);
  for (let i = 0; i < n; i++) {
    const o = {};
    for (let j = 0; j < tbl.columns.length; j++) o[tbl.columns[j]] = cols[j][i];
    out[i] = o;
  }
  return out;
}

const MISSING = `<span class="tag missing" title="Input absent — recorded as missing, never imputed">missing</span>`;
const tierBadge = (t) => `<span class="tier ${esc(t)}">Tier ${esc(t)}</span>`;
const layerTag = (l, label) =>
  `<span class="tag ${l}" title="Evidence layer ${l.toUpperCase()}">${esc(label || (l.toUpperCase() + " " + layerWord(l)))}</span>`;
function layerWord(l) {
  return { l1: "measured", l2: "derived", l3: "predicted", l4: "experimental" }[l] || "";
}

/* ---------------- charts (hand-rolled SVG) ---------------- */

function histogramSVG(values, { bins = 20, width = 560, height = 190, xLabel = "", color = "var(--l1)", from = 0, to = 1 } = {}) {
  if (!values.length) return `<p class="dim small">no values available</p>`;
  const counts = new Array(bins).fill(0);
  const span = to - from;
  for (const v of values) {
    if (v < from || v > to) continue;
    counts[Math.min(bins - 1, Math.floor(((v - from) / span) * bins))]++;
  }
  const pad = { l: 8, r: 8, t: 12, b: 26 };
  const max = Math.max(...counts);
  const iw = width - pad.l - pad.r, ih = height - pad.t - pad.b;
  const bw = iw / bins;
  const bars = [];
  for (let i = 0; i < bins; i++) {
    const h = max ? (counts[i] / max) * ih : 0;
    const lo = from + span * (i / bins);
    const hi = from + span * ((i + 1) / bins);
    bars.push(`<rect x="${(pad.l + i * bw).toFixed(1)}" y="${(pad.t + ih - h).toFixed(1)}" width="${(bw - 1.5).toFixed(1)}" height="${h.toFixed(1)}" fill="${color}" rx="1.5"><title>${counts[i]} records in [${lo.toFixed(3)}, ${hi.toFixed(3)})</title></rect>`);
  }
  const tick = (frac, txt) =>
    `<text class="axis" x="${(pad.l + frac * iw).toFixed(1)}" y="${height - 8}" text-anchor="middle">${esc(txt)}</text>`;
  return `<svg viewBox="0 0 ${width} ${height}" role="img" aria-label="histogram">${bars.join("")}
    ${tick(0, String(from))}${tick(0.5, (from + span / 2).toFixed(2))}${tick(1, String(to))}
    <text class="axis" x="${width - pad.r}" y="11" text-anchor="end">${esc(xLabel)} · peak bin ${nfmt(max)}</text>
    <line x1="${pad.l}" y1="${pad.t + ih}" x2="${width - pad.r}" y2="${pad.t + ih}" stroke="var(--line)"/></svg>`;
}

function barChartSVG(entries, { width = 560, labelW = 52, color = "var(--l2)" } = {}) {
  if (!entries.length) return "";
  const rowH = 22, padT = 4, height = padT + entries.length * rowH + 4;
  const iw = width - labelW - 60;
  const max = Math.max(...entries.map((e) => e.value)) || 1;
  const rows = entries.map((e, i) => {
    const w = (e.value / max) * iw;
    const y = padT + i * rowH;
    return `<text class="barlabel" x="${labelW - 6}" y="${y + 15}" text-anchor="end">${esc(e.label)}</text>
      <rect x="${labelW}" y="${y + 3}" width="${Math.max(0, w).toFixed(1)}" height="13" fill="${e.color || color}" rx="3"><title>${esc(e.label)}: ${nfmt(e.value)}</title></rect>
      <text class="barlabel" x="${(labelW + w + 6).toFixed(1)}" y="${y + 15}">${nfmt(e.value)}</text>`;
  });
  return `<svg viewBox="0 0 ${width} ${height}" role="img" aria-label="bar chart" style="max-width:${width}px">${rows.join("")}</svg>`;
}

function curveSVG(curve, { width = 560, height = 210, color = "var(--accent)", base = "var(--faint)" } = {}) {
  const k = curve.k_pct, found = curve.found_fraction, rnd = curve.random_expectation;
  const pad = { l: 38, r: 12, t: 12, b: 26 };
  const iw = width - pad.l - pad.r, ih = height - pad.t - pad.b;
  const to = (i, v) =>
    `${(pad.l + (i / (k.length - 1)) * iw).toFixed(1)},${(pad.t + ih - v * ih).toFixed(1)}`;
  const poly = (vals, col, dash) =>
    `<polyline fill="none" stroke="${col}" stroke-width="2" ${dash ? `stroke-dasharray="4 4"` : ""}
      points="${vals.map((v, i) => to(i, v)).join(" ")}"/>`;
  const dots = found.map((v, i) =>
    `<circle cx="${(pad.l + (i / (k.length - 1)) * iw).toFixed(1)}" cy="${(pad.t + ih - v * ih).toFixed(1)}" r="3" fill="${color}"><title>top ${k[i]}%: ${(100 * v).toFixed(1)}% of validated (random expect ${(100 * rnd[i]).toFixed(0)}%)</title></circle>`).join("");
  return `<svg viewBox="0 0 ${width} ${height}" style="max-width:${width}px" role="img" aria-label="rank at k curve">
    ${poly(rnd, base, true)}${poly(found, color, false)}${dots}
    ${k.map((kk, i) => `<text class="axis" x="${(pad.l + (i / (k.length - 1)) * iw).toFixed(1)}" y="${height - 8}" text-anchor="middle">${kk}%</text>`).join("")}
    <text class="axis" x="4" y="${pad.t + 4}">1.0</text>
    <line x1="${pad.l}" y1="${pad.t + ih}" x2="${width - pad.r}" y2="${pad.t + ih}" stroke="var(--line)"/></svg>`;
}

/* ---------------- state, data, routing ---------------- */

const CASES = [
  { id: "hcc1395", label: "HCC1395 / HCC1395BL (SEQC2 truth set)" },
  { id: "skcm_d9a3z1", label: "TCGA-D9-A3Z1 melanoma (open tier)" },
];

const state = {
  caseId: localStorage.getItem("nf.case") || "hcc1395",
  data: null,
  tesla: null,
  mu: { q: "", chrom: "", pepOnly: false, page: 1 },
  ca: { variant: "", allele: "", len: "", sort: "aff", tsp: false, page: 1 },
  rule: null,
  variant: null,
};

async function fetchJSON(url) {
  const res = await fetch(url);
  if (!res.ok) throw new Error(`${url} → HTTP ${res.status}`);
  if (url.endsWith(".json.gz")) {
    if (!("DecompressionStream" in window)) {
      throw new Error("this browser cannot decompress the exported projections " +
        "(DecompressionStream missing); use a current Chrome, Edge, Firefox or Safari");
    }
    const text = await new Response(
      new Response(res.arrayBuffer()).body.pipeThrough(new DecompressionStream("gzip"))
    ).text();
    return JSON.parse(text);
  }
  return res.json();
}

async function ensureCase(id) {
  if (state.caseId === id && state.data) return state.data;
  state.caseId = id;
  localStorage.setItem("nf.case", id);
  state.data = await fetchJSON(`data/case_${id}.json.gz`);
  state.variant = null;
  state.rule = null;
  state.mu = { q: "", chrom: "", pepOnly: false, page: 1 };
  state.ca = { variant: "", allele: "", len: "", sort: "aff", tsp: false, page: 1 };
  return state.data;
}
const ensureTesla = () => state.tesla || (state.tesla = fetchJSON("data/benchmark_tesla.json.gz"));

function parseHash() {
  const h = location.hash.replace(/^#\/?/, "");
  const [route, query] = h.split("?");
  return { route: route || "overview", params: new URLSearchParams(query || "") };
}

const VIEWS = [
  { id: "overview",    num: "",    label: "Overview",                render: viewOverview },
  { id: "fingerprint", num: "1",   label: "Tumour fingerprint",      render: viewFingerprint },
  { id: "mutations",   num: "2",   label: "Mutation explorer",       render: viewMutations },
  { id: "protein",     num: "3",   label: "Protein view",            render: viewProtein },
  { id: "candidates",  num: "4",   label: "Neoantigen candidates",   render: viewCandidates },
  { id: "prioritize",  num: "5",   label: "Prioritisation",          render: viewPrioritize },
  { id: "tesla",       num: "L4",  label: "TESLA benchmark",         render: viewTesla },
  { id: "workbench",   num: "WB",  label: "Analysis workbench",      render: viewWorkbench },
  { id: "report",      num: "6",   label: "Research report",         render: viewReport },
  { id: "construct",   num: "7",   label: "Construct (design)",      render: viewConstruct },
  { id: "rosie",       num: "",    label: "Rosie pathway mapping",   render: viewRosie },
];

function renderNav(active) {
  $("#navList").innerHTML = VIEWS.map((v) =>
    `<li><a href="#/${v.id}" class="${v.id === active ? "active" : ""}">
       ${v.num ? `<span class="navnum">${esc(v.num)}</span>` : ""}${esc(v.label)}</a></li>`).join("");
  $("#casePicker").value = state.caseId;
}

async function router() {
  const { route, params } = parseHash();
  const view = VIEWS.find((v) => v.id === route) || VIEWS[0];
  renderNav(view.id);
  app.innerHTML = `<div class="card"><p class="dim">loading exported artifact projections…</p></div>`;
  try {
    if (!state.tesla) { try { state.tesla = await ensureTesla(); } catch { state.tesla = null; } }
    await ensureCase(state.caseId);
    await view.render(params);
    window.scrollTo(0, 0);
  } catch (err) {
    app.innerHTML = `<div class="notice red"><span class="ttl">Data unavailable.</span> ${esc(err.message)}
      <p>This interface serves exported pipeline projections from <span class="mono">frontend/data/</span>.
      Regenerate them with <span class="mono">python -m neoforge.export_ui --work-dir … --cases-dir … --out frontend/data</span>
      and view over HTTP, e.g. <span class="mono">python -m http.server -d frontend</span>.</p></div>`;
  }
}
window.addEventListener("hashchange", router);

/* case picker */
$("#casePicker").innerHTML =
  CASES.map((c) => `<option value="${c.id}">${esc(c.label)}</option>`).join("");
$("#casePicker").addEventListener("change", async (e) => {
  await ensureCase(e.target.value).catch(() => {});
  router();
});

/* ---------------- shared case components ---------------- */

function statusTag(s) {
  return `<span class="tag ${s === "pass" ? "ok" : s === "warn" ? "warn" : "missing"}">${esc(s)}</span>`;
}

function caseHeader(d, title, question) {
  const b = d.case_bundle;
  const hla = d.ingest.hla;
  return `<h1>${esc(title)}</h1>
    <p class="sub">${esc(question)}</p>
    <p class="small"><span class="mono dim">${esc(b.display_name)}</span>
      ${layerTag("l1", "L1 measured")} ${layerTag("l2", "L2 derived")}
      ${d.present ? layerTag("l3", "L3 predicted") : ""}
      <span class="tag ${d.report ? "l4" : "missing"}">${esc(d.report ? d.report.layer : "chain ends at L2")}</span></p>
    ${hla && hla.status === "unavailable" ? `<div class="notice"><span class="ttl">Presentation stages mechanically refused for this case.</span><br>
      ${esc(hla.reason)} The chain therefore stops after protein-level derivation: L1 (measured) + L2 (derived) only —
      binding predictions must not be invented without the patient's HLA genotype.</div>` : ""}`;
}

function stageSourceLine(d, stage) {
  const rec = d.pipeline_stage_records[stage];
  return rec ? `<p class="small dim mono" title="${esc(rec.sha256)}">source_artifact: ${esc(rec.stage)}.${esc(rec.artifact_id)} · sha256 ${shortSha(rec.sha256)}</p>` : "";
}

function pagerHTML(page, pages, shown, total) {
  return `<div class="pager">
    <button data-page="prev" ${page <= 1 ? "disabled" : ""}>‹ prev</button>
    <span>page ${page} / ${pages} · showing ${nfmt(shown)} rows</span>
    <button data-page="next" ${page >= pages ? "disabled" : ""}>next ›</button></div>
    <p class="small dim">${nfmt(total)} rows match the current filters.</p>`;
}

/* ---------------- overview ---------------- */

function viewOverview() {
  const d = state.data;
  const b = d.case_bundle;
  const prof = d.report ? d.report.molecular_profile : {
    snv_truth_records: d.ingest.variant_counts.somatic_snv_positions,
    indel_truth_records: d.ingest.variant_counts.somatic_indel_positions,
    variants_with_peptides: d.translate.summary.variants_with_peptides,
    windows_derived: d.translate.summary.windows_derived,
  };
  const stages = ["ingest", "variants", "translate", "present", "prioritize", "report"];
  const sNames = { ingest: "Ingest", variants: "Variants", translate: "Peptides", present: "Present", prioritize: "Prioritize", report: "Report" };

  const t = state.tesla;
  const teslaCard = t ? `<div class="card">
      <h2 style="margin-top:0">Ranking-rule benchmark — TESLA ${layerTag("l4", "L4 experimental")}</h2>
      <p class="sub small">The candidate ranking rule was audited on the TESLA benchmark: ${nfmt(t.metrics.n_candidates)}
      pMHC predictions across ${t.metrics.n_patients} patients / ${t.metrics.n_alleles} alleles, compared against
      experimentally measured binding labels. Cohort-level evidence about the rule — never validation of any single case.</p>
      <div class="stats">
        <div class="stat"><div class="v">${num(t.metrics.ttif_top20pct, 3)}</div><div class="k">TTIF in top 20% · random expect ${num(t.metrics.ttif_random_expectation, 2)}</div></div>
        <div class="stat"><div class="v">${num(t.metrics.auprc, 3)}</div><div class="k">AUPRC · prevalence ${num(t.metrics.auprc_random_expectation, 3)}</div></div>
        <div class="stat"><div class="v">${nfmt(t.metrics.median_rank_of_validated)}</div><div class="k">median rank of the ${t.metrics.n_validated} validated</div></div>
        <div class="stat"><div class="v">${t.rank_rule_falsified ? "FALSIFIED" : "not falsified"}</div><div class="k">pre-registered verdict</div></div>
      </div>
      <p class="small"><a href="#/tesla">inspect the full benchmark page →</a></p></div>`
    : `<div class="notice"><span class="ttl">TESLA benchmark:</span> not exported (shown as pending until the harness has run — the interface never substitutes a number it cannot trace).</div>`;

  app.innerHTML = `
    <h1>Molecularly individual cancer, traced end to end</h1>
    <p class="sub">Every tumour carries its own set of mutations. NEOFORGE runs a real, deterministic computational
    pipeline on a real tumour–normal pair and lays the science out for inspection: what differs from healthy tissue,
    which altered protein fragments could be presented on this individual's HLA molecules, and why specific candidates
    are shortlisted for <em>further research</em>. It is inspired by the reported Rosie case, in which a computational
    analysis pathway fed expert-validated personalised vaccine research at UNSW — here we replicate that
    <em>computational research pathway</em>, never the treatment or its outcome.</p>

    <div class="card">
      <h2 style="margin-top:0">Research case — ${esc(b.display_name)}</h2>
      <p class="sub small">${esc(b.description)}</p>
      <div class="stats">
        <div class="stat"><div class="v">${nfmt(prof.snv_truth_records ?? 0)}</div><div class="k">somatic SNV records (L1)</div></div>
        <div class="stat"><div class="v">${nfmt(prof.indel_truth_records ?? 0)}</div><div class="k">somatic indel records (L1)</div></div>
        <div class="stat"><div class="v">${nfmt(prof.variants_with_peptides ?? 0)}</div><div class="k">variants with peptide windows (L2)</div></div>
        <div class="stat"><div class="v">${nfmt(prof.windows_derived ?? (d.translate ? d.translate.derived_rows_total : 0))}</div><div class="k">peptide windows derived (L2)</div></div>
        ${d.prioritize ? `<div class="stat"><div class="v">${nfmt(d.prioritize.candidates_prioritised)}</div><div class="k">candidates prioritised (L1+L3)</div></div>
        <div class="stat"><div class="v">${nfmt(d.report.shortlist.length)}</div><div class="k">report shortlist</div></div>` : ""}
      </div>
    </div>

    <div class="card">
      <h2 style="margin-top:0">The pipeline — one hash-verified chain${d.present ? "" : " (ends after peptide derivation: no same-patient HLA in the open tier)"}</h2>
      <div class="pipeline">
        ${stages.map((s, i) => {
          const rec = d.pipeline_stage_records[s];
          const off = !rec;
          const why = off ? (d.ingest.hla && d.ingest.hla.status === "unavailable" ? "refused — no HLA" : "not run") : esc(rec.artifact_id);
          return `<div class="pstage ${off ? "off" : ""}" title="${rec ? esc(rec.sha256) : "stage not run"}">
              <div class="n">STAGE ${i + 1}</div><div class="t">${esc(sNames[s])}</div>
              <div class="h">${why}</div></div>
            ${i < stages.length - 1 ? `<span class="parrow">→</span>` : ""}`;
        }).join("")}
      </div>
      <p class="small"><a href="#/workbench">inspect artifacts + SHA-256 provenance in the workbench →</a></p>
    </div>

    <div class="two-col">
      <div class="card">
        <h2 style="margin-top:0">The demo path — one continuous research journey</h2>
        <ul>
          <li>Step 1 — <a href="#/fingerprint">tumour fingerprint</a>: the measured molecular differences vs the matched normal.</li>
          <li>Step 2 — <a href="#/mutations">mutation explorer</a>: individual variants and their evidence classes.</li>
          <li>Step 3 — <a href="#/protein">protein view</a>: codon-exact reference vs mutant peptide windows.</li>
          <li>Step 4 — <a href="#/candidates">neoantigen candidates</a>: predicted HLA binding, explicitly tagged L3.</li>
          <li>Step 5 — <a href="#/prioritize">prioritisation</a>: a fully inspectable rule with a live playground.</li>
          <li>L4 — <a href="#/tesla">TESLA benchmark</a>: external, cohort-level audit of the ranking rule.</li>
          <li>WB — <a href="#/workbench">workbench</a> · 6 — <a href="#/report">research report</a>: provenance, uncertainty, export.</li>
        </ul>
          <li>Step 7 — <a href="#/construct">construct (design)</a>: an in-silico polyepitope research design — every epitope choice ledgered, never a manufacturable claim.</li>
          <li><a href="#/rosie">Rosie pathway mapping</a>: the reported dog case, step by step, mapped to what NEOFORGE does and does not do.</li>
      </div>
      <div class="notice" style="margin:14px 0">
        <span class="ttl">Deliberately absent from every screen</span>
        <ul class="small" style="margin:6px 0 0;padding-left:18px">
          <li>any claim of efficacy, diagnosis, or treatment selection</li>
          <li>any immunogenicity prediction</li>
          <li>cross-patient data blending — every input is same-patient</li>
          <li>composite biological scores or invented confidences</li>
        </ul>
      </div>
    </div>

    ${teslaCard}`;
}

/* ---------------- step 1: fingerprint ---------------- */

function viewFingerprint() {
  const d = state.data;
  const ing = d.ingest;
  const snv = ing.snv_table;
  const tvafVals = (snv.tvaf || []).filter((v) => v !== null);
  const nvafVals = (snv.nvaf || []).filter((v) => v !== null);
  const median = (a) => { const s = [...a].sort((x, y) => x - y); return s.length ? s[Math.floor(s.length / 2)] : null; };
  const chromOrder = (c) => {
    const base = String(c).replace(/^chr/, "");
    const m = { X: 23, Y: 24, M: 25, MT: 25 };
    const n = parseInt(base, 10);
    return Number.isFinite(n) ? n : (m[base.toUpperCase()] ?? 99);
  };
  const chromCounts = {};
  for (const c of snv.chrom) { const base = "chr" + String(c).replace(/^chr/, ""); chromCounts[base] = (chromCounts[base] || 0) + 1; }
  const chromEntries = Object.entries(chromCounts)
    .sort((a, b) => chromOrder(a[0]) - chromOrder(b[0]))
    .map(([label, value]) => ({ label: label.replace("chr", ""), value }));

  const vafSource = ing.maf
    ? "TVAF computed from MAF allele depths: tvaf = t_alt_count / t_depth (the single documented transform); anomalies counted, not dropped"
    : "TVAF read from the VCF INFO/TVAF field (SEQC2 truth-set annotation)";

  app.innerHTML = `
    ${caseHeader(d, "Tumour fingerprint", "What makes this tumour molecularly different from its matched normal — measured, not predicted?")}
    ${stageSourceLine(d, "ingest")}

    <div class="stats">
      <div class="stat"><div class="v">${nfmt(ing.variant_counts.somatic_snv_positions)}</div><div class="k">high-confidence somatic SNVs</div></div>
      <div class="stat"><div class="v">${nfmt(ing.variant_counts.somatic_indel_positions)}</div><div class="k">somatic indels</div></div>
      ${d.report ? `<div class="stat"><div class="v">${nfmt(d.report.molecular_profile.transcripts_hit)}</div><div class="k">transcripts hit (L2)</div></div>
      <div class="stat"><div class="v">${nfmt(d.report.molecular_profile.peptide_windows)}</div><div class="k">peptide windows derived (L2)</div></div>` : ""}
    </div>

    <div class="two-col">
      <div class="card">
        <h3 style="margin-top:0">Tumour variant allele fraction ${layerTag("l1", "L1 measured")}</h3>
        <p class="small dim">Distribution of ${nfmt(tvafVals.length)} SNVs with computable TVAF
        (median ${num(median(tvafVals), 3)}). High values mean the mutant allele dominates the tumour tissue; the
        matched normal sits near zero everywhere (median ${num(median(nvafVals), 4)}) — that separation is what makes
        these records tumour-specific observations rather than assumptions.</p>
        <div class="chart-wrap">${histogramSVG(tvafVals, { xLabel: "tumour VAF", color: "var(--l1)", bins: 20 })}</div>
        <p class="small dim mono">${esc(vafSource)}</p>
      </div>
      <div class="card">
        <h3 style="margin-top:0">Mutations per chromosome</h3>
        <div class="chart-wrap">${barChartSVG(chromEntries)}</div>
      </div>
    </div>

    <div class="card">
      <h3 style="margin-top:0">Input validation checks</h3>
      <table><tbody>
        ${ing.checks.map((c) => `<tr><td style="width:56px">${statusTag(c.status)}</td>
          <td class="mono" style="width:340px">${esc(c.id)}</td><td class="wrap">${esc(c.detail)}</td></tr>`).join("")}
      </tbody></table>
    </div>

    <div class="card">
      <h3 style="margin-top:0">Samples, HLA type & expression — the inputs the rest of the chain stands on</h3>
      <dl class="kv">
        ${ing.maf ? `<dt>tumour sample</dt><dd>${esc(ing.maf.tumour_sample)}</dd>
                    <dt>matched normal</dt><dd>${esc(ing.maf.normal_sample)}</dd>` : ""}
        ${ing.hla.allele_statuses ? `<dt>HLA class I ${layerTag("l1", "measured typing")}</dt>
                    <dd>${ing.hla.allele_statuses.map((a) => esc(a.normalised)).filter(Boolean).join(", ")}</dd>` : `<dt>HLA class I</dt><dd>${MISSING}</dd>`}
        ${ing.expression ? `<dt>expression</dt><dd class="wrap">${esc(ing.expression.matrix || (ing.expression.source && ing.expression.source.path) || "")} · ${nfmt(ing.expression.n_genes)} genes</dd>` : ""}
      </dl>
      <p class="small"><a href="#/workbench">full checksums, versions, licences → workbench</a></p>
    </div>`;
}

/* ---------------- step 2: mutation explorer ---------------- */

function mutationRows(d) {
  const sel = d.variants ? rowsOf(d.variants.selected) : [];
  const winCounts = {};
  if (d.translate) for (const r of rowsOf(d.translate.derived)) winCounts[r.variant_key] = (winCounts[r.variant_key] || 0) + 1;
  const cand = {};
  if (d.prioritize) for (const c of rowsOf(d.prioritize.candidates)) cand[c.variant_key] = c;
  const tierRank = { A: 0, B: 1, C: 2 };
  const out = sel.map((v) => {
    const c = cand[v.variant_key];
    return { ...v, windows: winCounts[v.variant_key] || 0, tier: c ? c.tier : null, cand: c };
  });
  out.sort((a, b) =>
    ((b.windows > 0) - (a.windows > 0)) ||
    ((a.tier ? tierRank[a.tier] : 3) - (b.tier ? tierRank[b.tier] : 3)) ||
    (a.cand && b.cand ? a.cand.mhcflurry_affinity_nm - b.cand.mhcflurry_affinity_nm : 0));
  return out;
}

function evidenceTags(r, d) {
  const tags = [layerTag("l1", "measured VAF/depth")];
  if (r.windows > 0) tags.push(layerTag("l2", "peptides derived"));
  if (d.present && r.cand) tags.push(layerTag("l3", "binding predicted"));
  return tags.join(" ");
}

function viewMutations() {
  const d = state.data;
  const rows = mutationRows(d);
  const chroms = [...new Set(rows.map((r) => r.chrom))].sort();
  const tiers = d.prioritize ? ["A", "B", "C"] : [];
  const mu = state.mu;

  app.innerHTML = `
    ${caseHeader(d, "Mutation explorer", "Which somatic variants exist, what evidence accompanies each, and which of them actually gave rise to altered protein fragments?")}
    ${stageSourceLine(d, "variants")}

    <div class="controls">
      <input id="mu_q" placeholder="filter variant… e.g. chr7:551" value="${esc(mu.q)}">
      <select id="mu_chrom"><option value="">all chromosomes</option>
        ${chroms.map((c) => `<option ${mu.chrom === c ? "selected" : ""} value="${esc(c)}">${esc(c)}</option>`).join("")}</select>
      ${d.prioritize ? `<select id="mu_tier"><option value="">all tiers</option>
        ${tiers.map((t) => `<option ${mu.tier === t ? "selected" : ""} value="${t}">Tier ${t}</option>`).join("")}</select>` : ""}
      <label><input type="checkbox" id="mu_pep" ${mu.pepOnly ? "checked" : ""}> with derived peptides only</label>
    </div>
    <div id="mu_table"></div>
    <div id="mu_drawer"></div>`;

  const renderTable = () => {
    const q = ($("#mu_q").value || "").trim();
    const chrom = $("#mu_chrom").value;
    const tier = $("#mu_tier") ? $("#mu_tier").value : "";
    const pepOnly = $("#mu_pep").checked;
    Object.assign(mu, { q, chrom, tier, pepOnly });

    let filtered = rows;
    if (q) filtered = filtered.filter((r) => r.variant_key.includes(q));
    if (chrom) filtered = filtered.filter((r) => r.chrom === chrom);
    if (tier) filtered = filtered.filter((r) => r.tier === tier);
    if (pepOnly) filtered = filtered.filter((r) => r.windows > 0);

    const PAGE = 40;
    const pages = Math.max(1, Math.ceil(filtered.length / PAGE));
    mu.page = Math.min(Math.max(1, mu.page), pages);
    const slice = filtered.slice((mu.page - 1) * PAGE, mu.page * PAGE);

    $("#mu_table").innerHTML = `
      <table><thead><tr><th>variant</th><th>tumour VAF (L1)</th><th>normal VAF (L1)</th>
        <th class="num">depth</th>${d.prioritize ? "<th>tier</th>" : ""}<th>peptide windows (L2)</th><th>evidence</th></tr></thead>
      <tbody>${slice.map((r) => `
        <tr class="rowlink" data-key="${esc(r.variant_key)}">
          <td>${esc(r.variant_key)}</td>
          <td>${r.tvaf === null || r.tvaf === undefined ? MISSING : num(r.tvaf, 3)}</td>
          <td>${r.nvaf === null || r.nvaf === undefined ? MISSING : num(r.nvaf, 3)}</td>
          <td class="num">${r.tumor_depth === null || r.tumor_depth === undefined ? "—" : nfmt(r.tumor_depth)}</td>
          ${d.prioritize ? `<td>${r.tier ? tierBadge(r.tier) : MISSING}</td>` : ""}
          <td>${r.windows ? nfmt(r.windows) : `<span class="dim">0</span>`}</td>
          <td>${evidenceTags(r, d)}</td>
        </tr>`).join("")}</tbody></table>
      ${pagerHTML(mu.page, pages, slice.length, filtered.length)}`;

    document.querySelectorAll("#mu_table .rowlink").forEach((tr) =>
      tr.addEventListener("click", () => openMutationDrawer(tr.dataset.key, d)));
    document.querySelectorAll("#mu_table .pager button").forEach((b) =>
      b.addEventListener("click", () => { mu.page += b.dataset.page === "next" ? 1 : -1; renderTable(); }));
  };
  [["#mu_q", "input"], ["#mu_chrom", "change"], ["#mu_pep", "change"], ["#mu_tier", "change"]]
    .forEach(([id, ev]) => { const el = $(id); if (el) { el.addEventListener(ev, () => { mu.page = 1; renderTable(); }); if (id === "#mu_tier" || id === "#mu_pep" || id === "#mu_chrom") el.addEventListener(ev, renderTable); } });
  renderTable();
}

function openMutationDrawer(key, d) {
  const r = rowsOf(d.variants.selected).find((x) => x.variant_key === key);
  const c = d.prioritize ? rowsOf(d.prioritize.candidates).find((x) => x.variant_key === key) : null;
  const wins = d.translate ? rowsOf(d.translate.derived).filter((x) => x.variant_key === key) : [];
  const onShortlist = d.report ? d.report.shortlist.some((s) => s.variant_key === key) : false;

  $("#mu_drawer").innerHTML = `
    <div class="detail">
      <button class="close" aria-label="close">✕</button>
      <h2 style="margin-top:0" class="mono" style="font-size:17px">${esc(key)}</h2>
      ${onShortlist ? `<div class="notice good small" style="margin:8px 0"><span class="ttl">On the research-report shortlist.</span></div>` : ""}
      <h3>Measured (L1)</h3>
      <dl class="kv">
        <dt>tumour VAF</dt><dd>${r.tvaf === null || r.tvaf === undefined ? MISSING : num(r.tvaf, 3)}</dd>
        <dt>normal VAF</dt><dd>${r.nvaf === null || r.nvaf === undefined ? MISSING : num(r.nvaf, 3)}</dd>
        <dt>tumour depth</dt><dd>${r.tumor_depth === null || r.tumor_depth === undefined ? "—" : nfmt(r.tumor_depth)}</dd>
      </dl>
      ${c ? `<h3>Prioritisation (L1 + L3 rule)</h3>
      <dl class="kv">
        <dt>tier</dt><dd>${tierBadge(c.tier)}</dd>
        <dt>best peptide</dt><dd>${esc(c.best_peptide)} · ${esc(c.best_allele)}</dd>
        <dt>pred. affinity (L3)</dt><dd>${num(c.mhcflurry_affinity_nm, 1)} nM (reference ${num(c.reference_mhcflurry_affinity_nm, 1)} nM)</dd>
        <dt>presentation score (L3)</dt><dd>${num(c.mhcflurry_presentation_score, 3)}</dd>
        <dt>expression TPM (L1)</dt><dd>${c.expression_tpm === null || c.expression_tpm === undefined ? MISSING : num(c.expression_tpm, 2)}</dd>
        <dt>gate flags</dt><dd>${c.gate_flags.length ? esc(c.gate_flags.join(", ")) : "—"}</dd>
      </dl>` : ""}
      <h3>Derived protein windows (L2)</h3>
      ${wins.length ? wins.slice(0, 6).map((w) => `
        <p class="small mono" style="margin:6px 0 2px">${esc(w.transcript_id)} · p.${w.protein_pos} ${esc(w.ref_aa)}→${esc(w.alt_aa)} · ${w.peptide_length}-mer</p>
        ${seqTracksHTML(w)}`).join("")
        : `<p class="small dim">no protein-coding peptide windows derived for this variant (documented in the translate artifact — a classification, not an omission)</p>`}
      <div style="margin-top:14px;display:flex;gap:8px;flex-wrap:wrap">
        ${wins.length ? `<a href="#/protein?variant=${encodeURIComponent(key)}"><button class="controls">protein view →</button></a>` : ""}
        ${d.present ? `<a href="#/candidates?variant=${encodeURIComponent(key)}"><button class="controls">predictions →</button></a>` : ""}
      </div>
    </div>`;
  const close = () => { $("#mu_drawer").innerHTML = ""; document.removeEventListener("keydown", onKey); };
  const onKey = (e) => { if (e.key === "Escape") close(); };
  document.addEventListener("keydown", onKey);
  $(".detail .close").addEventListener("click", close);
}

/** reference vs mutant tracks with the altered residues highlighted */
function seqTracksHTML(w) {
  const ref = w.ref_window_seq, mut = w.window_seq;
  const hlRef = [...ref].map((ch, i) => (mut[i] !== undefined && ch !== mut[i]) ? `<b>${esc(ch)}</b>` : esc(ch)).join("");
  const hlMut = [...ref].map((ch, i) => (mut[i] !== undefined && ch !== mut[i]) ? `<b>${esc(mut[i])}</b>` : esc(mut[i] ?? "")).join("");
  const multi = [...ref].some((ch, i) => mut[i] !== undefined && ch !== mut[i]);
  return `<div class="seqtrack na"><span class="winlabel">ref&nbsp;&nbsp;&nbsp;</span> ${hlRef}</div>
     <div class="seqtrack"><span class="winlabel">mut&nbsp;&nbsp;&nbsp;</span> ${hlMut}</div>
     <div class="winlabel">${w.peptide_length}-mer · window starts at protein position ${w.window_start + 1}${multi ? " · altered residues highlighted" : " · identical to reference at peptide ends"}</div>`;
}

/* ---------------- step 3: protein view ---------------- */

function viewProtein(params) {
  const d = state.data;
  const key = params.get("variant") || state.variant;
  const derived = d.translate ? rowsOf(d.translate.derived) : [];

  if (!key || !derived.some((w) => w.variant_key === key)) {
    const shortlist = d.report ? d.report.shortlist.map((s) => s.variant_key) : [];
    const withWins = [...new Set(derived.map((w) => w.variant_key))].sort();
    const keys = shortlist.filter((k) => withWins.includes(k));
    app.innerHTML = `
      ${caseHeader(d, "Protein view", "How does one nucleotide change become an altered protein fragment?")}
      <div class="notice blue"><span class="ttl">Pick a variant.</span> Everything below is a codon-exact sequence
      derivation — a conceptual diagram, not a predicted 3-D structure.</div>
      <div class="card">
        <h3 style="margin-top:0">Shortlisted variants with derived windows</h3>
        <div style="display:flex;flex-direction:column;gap:4px;max-width:640px">
        ${keys.map((k) => `<a class="mono" href="#/protein?variant=${encodeURIComponent(k)}" style="text-decoration:none;padding:6px 10px;border:1px solid var(--line-soft);border-radius:8px">${esc(k)}</a>`).join("")}
        </div>
        <p class="small dim" style="margin-top:10px">${nfmt(withWins.length)} variants have derived windows in total
        (${esc(d.translate ? d.translate.export_scope : "no translate artifact")});
        the ${keys.length} shown are the report shortlist. Any of them can be opened by URL.</p>
      </div>`;
    return;
  }

  state.variant = key;
  const wins = derived.filter((w) => w.variant_key === key);
  const byTranscript = new Map();
  for (const w of wins) {
    if (!byTranscript.has(w.transcript_id)) byTranscript.set(w.transcript_id, []);
    byTranscript.get(w.transcript_id).push(w);
  }
  for (const list of byTranscript.values()) list.sort((a, b) => a.peptide_length - b.peptide_length || a.window_start - b.window_start);

  const noHLA = d.ingest.hla && d.ingest.hla.status === "unavailable";
  app.innerHTML = `
    ${caseHeader(d, "Protein view", "Codon-exact reference vs mutant peptide windows, derived from the annotated transcript for one variant.")}
    ${stageSourceLine(d, "translate")}
    <div class="notice blue small"><span class="ttl">Conceptual diagram.</span> One-letter amino-acid sequence windows
    derived codon-exactly from the reference genome annotation (L2). Not predicted 3-D structures; no structural claim.</div>

    <div class="card">
      <h2 style="margin-top:0" class="mono" style="font-size:18px">${esc(key)}</h2>
      <p class="small dim">${nfmt(wins.length)} windows across ${byTranscript.size} transcript(s) ·
        altered residues highlighted in the read-out below</p>
      ${[...byTranscript.entries()].map(([tid, list]) => `
        <h3>${esc(tid)} · p.${list[0].protein_pos} ${esc(list[0].ref_aa)}→${esc(list[0].alt_aa)} · strand ${esc(list[0].strand)}</h3>
        ${list.map((w) => `
          <div style="margin:10px 0;padding:10px 14px;background:var(--panel-2);border:1px solid var(--line-soft);border-radius:10px">
            ${seqTracksHTML(w)}
          </div>`).join("")}`).join("")}
    </div>
    <p class="small">${d.present
      ? `<a href="#/candidates?variant=${encodeURIComponent(key)}">see L3 binding predictions for these windows →</a>`
      : (noHLA ? `<span class="dim">No binding predictions exist for this case — the patient's HLA genotype is not
          available in the open tier (declared missing, not imputed), so the chain stops here mechanically.</span>` : "")}</p>`;
}

/* ---------------- step 4: candidates ---------------- */

function viewCandidates(params) {
  const d = state.data;
  if (!d.present) {
    app.innerHTML = caseHeader(d, "Neoantigen candidates", "") +
      `<div class="notice"><span class="ttl">Stage not run — by design.</span> This case has no same-patient HLA
      genotype in the open data tier, so the pipeline mechanically refuses to invent binding predictions.
      Every prior stage (L1 measured, L2 derived) remains fully inspectable in the
      <a href="#/mutations">mutation explorer</a>.</div>`;
    return;
  }
  const preds = d.present.predictions;
  const rows = rowsOf(preds);
  const keys = [...new Set(rows.map((r) => r.variant_key))].sort();
  const alleles = [...new Set(rows.map((r) => r.allele))].sort();
  const lengths = [...new Set(rows.map((r) => r.peptide_length))].sort((a, b) => a - b);
  const ca = state.ca;
  if (params.get("variant") && keys.includes(params.get("variant"))) ca.variant = params.get("variant");

  app.innerHTML = `
    ${caseHeader(d, "Neoantigen candidates", "For shortlisted variants: every peptide window × HLA allele × length, with predicted binding — explicitly the prediction layer, never measured immunogenicity.")}
    ${stageSourceLine(d, "present")}
    <p class="small dim">${nfmt(preds.rows_exported)} prediction rows exported (scope: ${esc(preds.export_scope)}) of
      ${nfmt(preds.rows_total)} computed ${layerTag("l3", "MHCflurry predictions")} — the full table stays hash-pinned in the artifact store and listed in the workbench.</p>

    <div class="controls">
      <select id="ca_var"><option value="">all exported variants</option>
        ${keys.map((k) => `<option ${ca.variant === k ? "selected" : ""} value="${esc(k)}">${esc(k)}</option>`).join("")}</select>
      <select id="ca_allele"><option value="">all alleles</option>
        ${alleles.map((a) => `<option value="${esc(a)}">${esc(a)}</option>`).join("")}</select>
      <select id="ca_len"><option value="">all lengths</option>
        ${lengths.map((l) => `<option value="${l}">${l}-mer</option>`).join("")}</select>
      <select id="ca_sort">
        <option value="aff" ${ca.sort === "aff" ? "selected" : ""}>sort: predicted affinity (strongest first)</option>
        <option value="pct" ${ca.sort === "pct" ? "selected" : ""}>sort: presentation percentile (best first)</option>
        <option value="delta" ${ca.sort === "delta" ? "selected" : ""}>sort: Δ vs reference</option></select>
      <label><input type="checkbox" id="ca_tsp" ${ca.tsp ? "checked" : ""}> mutant predicted to bind better than reference</label>
    </div>
    <div id="ca_table"></div>

    <div class="card"><h3 style="margin-top:0">Reading the columns</h3>
      <p class="small dim"><span class="mono">affinity</span> — predicted MHC binding in nM (lower = stronger) ${layerTag("l3")}.
      <span class="mono">presentation score / percentile</span> — MHCflurry's antigen-processing-and-presentation model:
      a different biological quantity from binding affinity, kept in a separate column for exactly that reason.
      <span class="mono">reference nM</span> — the paired healthy homologous peptide: for personalised research a
      candidate is only meaningful if it plausibly binds <em>better</em> than its normal counterpart.
      None of these quantities involve T cells, and nothing on this page predicts an immune response.</p>
      <p class="small dim">${esc(d.present.parameters.predictor)} (MHCflurry ${esc(d.present.parameters.mhcflurry_version)}) ·
      patient alleles: ${esc((d.present.parameters.alleles || []).join(", "))}</p></div>`;

  const renderTable = () => {
    ca.allele = $("#ca_allele").value; ca.len = $("#ca_len").value;
    ca.sort = $("#ca_sort").value; ca.tsp = $("#ca_tsp").checked; ca.variant = $("#ca_var").value;
    let rs = rows;
    if (ca.variant) rs = rs.filter((r) => r.variant_key === ca.variant);
    if (ca.allele) rs = rs.filter((r) => r.allele === ca.allele);
    if (ca.len) rs = rs.filter((r) => r.peptide_length === Number(ca.len));
    if (ca.tsp) rs = rs.filter((r) => r.reference_mhcflurry_affinity_nm !== null &&
      r.mhcflurry_affinity_nm < r.reference_mhcflurry_affinity_nm && r.peptide !== r.reference_peptide);
    rs = [...rs].sort((a, b) =>
      ca.sort === "aff" ? a.mhcflurry_affinity_nm - b.mhcflurry_affinity_nm :
      ca.sort === "pct" ? a.mhcflurry_presentation_percentile - b.mhcflurry_presentation_percentile :
      (a.mhcflurry_affinity_nm - a.reference_mhcflurry_affinity_nm) - (b.mhcflurry_affinity_nm - b.reference_mhcflurry_affinity_nm));

    const PAGE = 40;
    const pages = Math.max(1, Math.ceil(rs.length / PAGE));
    ca.page = Math.min(Math.max(1, ca.page), pages);
    const slice = rs.slice((ca.page - 1) * PAGE, ca.page * PAGE);

    $("#ca_table").innerHTML = `
      <table><thead><tr><th>variant</th><th>peptide</th><th>allele</th><th>len</th>
        <th class="num">affinity nM (L3)</th><th class="num">presentation score</th><th class="num">percentile</th>
        <th class="num">reference nM</th><th>tumour-specific?</th></tr></thead>
      <tbody>${slice.map((r) => {
        const tighter = r.reference_mhcflurry_affinity_nm !== null && r.mhcflurry_affinity_nm < r.reference_mhcflurry_affinity_nm;
        const delta = tighter || r.reference_mhcflurry_affinity_nm !== null
          ? r.mhcflurry_affinity_nm - r.reference_mhcflurry_affinity_nm : null;
        return `<tr><td class="mono small">${esc(r.variant_key)}</td><td><b>${esc(r.peptide)}</b></td>
          <td class="small">${esc(r.allele)}</td><td class="num">${r.peptide_length}</td>
          <td class="num">${num(r.mhcflurry_affinity_nm, 1)}</td>
          <td class="num">${num(r.mhcflurry_presentation_score, 3)}</td>
          <td class="num">${num(r.mhcflurry_presentation_percentile, 2)}</td>
          <td class="num">${r.reference_mhcflurry_affinity_nm === null || r.reference_mhcflurry_affinity_nm === undefined ? MISSING : num(r.reference_mhcflurry_affinity_nm, 1)}</td>
          <td>${r.peptide === r.reference_peptide
            ? `<span class="dim small">identical — not tumour-specific</span>`
            : (delta === null ? MISSING
               : tighter ? `<span style="color:var(--l1)">✓ Δ ${num(delta, 1)} nM</span>`
                         : `<span class="dim small">Δ +${num(delta, 1)} nM</span>`)}</td></tr>`;
      }).join("")}</tbody></table>
      ${pagerHTML(ca.page, pages, slice.length, rs.length)}`;
    document.querySelectorAll("#ca_table .pager button").forEach((b) =>
      b.addEventListener("click", () => { ca.page += b.dataset.page === "next" ? 1 : -1; renderTable(); }));
  };
  [["#ca_var", "change"], ["#ca_allele", "change"], ["#ca_len", "change"], ["#ca_sort", "change"], ["#ca_tsp", "change"]]
    .forEach(([id, ev]) => $(id).addEventListener(ev, () => { ca.page = 1; renderTable(); }));
  renderTable();
}

/* ---------------- step 5: prioritisation + playground ---------------- */

/** Exact client-side reproduction of stages/prioritize.py's tier rule,
    including its missing-gate semantics (null never passes as satisfied). */
function recomputeTier(c, rule) {
  const strong = c.mhcflurry_affinity_nm < rule.nm;
  const expressed = c.expression_tpm === null || c.expression_tpm === undefined ? null : c.expression_tpm >= rule.tpm;
  const passesVaf = c.tvaf === null || c.tvaf === undefined ? null : c.tvaf >= rule.vaf;
  const tier = (strong && expressed === true && passesVaf !== false) ? "A" : (strong ? "B" : "C");
  return { strong, expressed, passesVaf, tier };
}

function viewPrioritize() {
  const d = state.data;
  if (!d.prioritize) {
    app.innerHTML = caseHeader(d, "Prioritisation", "") +
      `<div class="notice"><span class="ttl">Stage not run — by design.</span> Prioritisation joins L3 binding
      predictions with measured evidence; no HLA genotype means no predictions means no tiers.
      See the <a href="#/mutations">mutation explorer</a> for this case's L1/L2 evidence.</div>`;
    return;
  }
  const p = d.prioritize;
  const pr = p.parameters;
  const cands = rowsOf(p.candidates);
  state.rule = state.rule || { nm: pr.strong_binder_nm, tpm: pr.expressed_tpm, vaf: pr.min_tumour_vaf };
  const rule = state.rule;

  /* baseline self-check: the published rule must reproduce the artifact tiers */
  const countsOf = (assigned) => {
    const c = { A: 0, B: 0, C: 0 };
    for (const r of assigned.values()) c[r.tier]++;
    return c;
  };
  const baselineAssigned = new Map(cands.map((c) => [c.variant_key, recomputeTier(c, rule)]));
  const baselineCounts = countsOf(baselineAssigned);
  const reproduces = JSON.stringify(baselineCounts) === JSON.stringify({ A: p.tier_counts.A || 0, B: p.tier_counts.B || 0, C: p.tier_counts.C || 0 });

  const gateExplain = (r) => [
    r.strong ? `affinity < ${rule.nm} nM` : `affinity ≥ ${rule.nm} nM (not a strong binder under your threshold)`,
    r.expressed === null ? "TPM missing → gate recorded missing" : (r.expressed ? `TPM ≥ ${rule.tpm}` : `TPM < ${rule.tpm}`),
    r.passesVaf === null ? "VAF missing → gate recorded missing" : (r.passesVaf ? `VAF ≥ ${rule.vaf}` : `VAF < ${rule.vaf}`),
  ].join("; ");

  const tc = p.tier_counts;
  const totalC = cands.length;
  app.innerHTML = `
    ${caseHeader(d, "Prioritisation", "How does the pipeline rank candidates — explicitly, with documented thresholds and no invented composite scores?")}
    ${stageSourceLine(d, "prioritize")}

    <div class="card">
      <h3 style="margin-top:0">The published rule — no weighted score anywhere in the pipeline</h3>
      <div class="rulebox">strong binder      best predicted affinity  <  ${pr.strong_binder_nm} nM    (L3)
expressed          measured TPM             >=  ${pr.expressed_tpm}             (L1) — gate recorded missing when the source lacks the gene
passes VAF gate    measured tumour VAF      >=  ${pr.min_tumour_vaf}            (L1) — gate recorded missing when depth never allowed TVAF

Tier A = strong binder AND expressed AND VAF-gate not failed
Tier B = strong binder (expression or VAF missing / sub-threshold)
Tier C = everything else
missing gate input is never silently treated as passed
ordering: tier, then ascending best predicted affinity</div>
      <p class="small dim">threshold provenance: ${esc((pr.threshold_sources || []).join(" · "))} ·
      expression source: ${esc(pr.expression_source)}</p>
    </div>

    <div class="card">
      <h3 style="margin-top:0">Tier outcome — ${nfmt(totalC)} candidates</h3>
      ${barChartSVG(["A", "B", "C"].map((t) => ({ label: "Tier " + t, value: tc[t] || 0,
        color: t === "A" ? "var(--l1)" : t === "B" ? "var(--l3)" : "#2c3a4c" })), { width: 420 })}
      <p class="small dim mono">client-side recomputation at the published defaults reproduces the artifact tiers:
      ${reproduces
        ? `<span style="color:var(--l1)">✓ yes (${tc.A || 0}/${tc.B || 0}/${tc.C || 0})</span>`
        : `<span style="color:var(--bad)">✗ NO (${baselineCounts.A}/${baselineCounts.B}/${baselineCounts.C}) — do not trust the playground for this export</span>`}</p>
    </div>

    <div class="card">
      <h3 style="margin-top:0">Rule playground — move any threshold, watch membership re-derive</h3>
      <p class="small dim">The same fixed rule, evaluated in your browser over the exported candidate features,
      preserving the missing-gate semantics. The artifact itself is immutable: the playground never writes back,
      and the counts above keep showing the artifact values for comparison.</p>
      <div class="playgrid">
        <div class="playcard"><label>strong-binder cut</label>
          <input type="range" id="pl_nm" min="50" max="2000" step="50" value="${rule.nm}">
          <div class="cur">affinity &lt; <span id="pl_nm_v">${rule.nm}</span> nM</div></div>
        <div class="playcard"><label>expression gate</label>
          <input type="range" id="pl_tpm" min="0" max="10" step="0.5" value="${rule.tpm}">
          <div class="cur">TPM ≥ <span id="pl_tpm_v">${rule.tpm}</span></div></div>
        <div class="playcard"><label>tumour-VAF gate</label>
          <input type="range" id="pl_vaf" min="0" max="0.9" step="0.05" value="${rule.vaf}">
          <div class="cur">VAF ≥ <span id="pl_vaf_v">${rule.vaf}</span></div></div>
      </div>
      <h3>Tier counts under your rule <span class="dim small">(artifact values in parentheses)</span></h3>
      <div id="pl_out_counts"></div>
      <h3>Candidates changing tier</h3>
      <div id="pl_out_diff"></div>
    </div>

    <div class="card">
      <h3 style="margin-top:0">Shortlist carried into the research report</h3>
      <table><thead><tr><th>variant</th><th>tier</th><th>best peptide</th><th>allele</th>
        <th class="num">affinity nM (L3)</th><th class="num">TVAF (L1)</th><th class="num">TPM (L1)</th></tr></thead><tbody>
        ${d.report.shortlist.map((s) => `<tr><td class="mono small">${esc(s.variant_key)}</td><td>${tierBadge(s.tier)}</td>
          <td><b>${esc(s.best_peptide)}</b></td><td class="small">${esc(s.best_allele)}</td>
          <td class="num">${num(s.mhcflurry_affinity_nm, 1)}</td>
          <td class="num">${num(s.tvaf, 3)}</td><td class="num">${num(s.expression_tpm, 1)}</td></tr>`).join("")}
      </tbody></table>
      <p class="small"><a href="#/report">open the full research report →</a></p>
    </div>`;

  const applyPlay = () => {
    const r = { nm: Number($("#pl_nm").value), tpm: Number($("#pl_tpm").value), vaf: Number($("#pl_vaf").value) };
    state.rule = r;
    const assigned = new Map(cands.map((c) => [c.variant_key, recomputeTier(c, r)]));
    const c = countsOf(assigned);
    $("#pl_out_counts").innerHTML = `
      <table><tbody>
        ${["A", "B", "C"].map((t) => `<tr><td style="width:90px">${tierBadge(t)}</td>
          <td class="num" style="width:80px">${nfmt(c[t])}</td>
          <td class="num dim small" style="width:90px">(artifact: ${nfmt(tc[t] || 0)})</td>
          <td>${c[t] === (tc[t] || 0) ? `<span class="tag ok">unchanged</span>` : `<span class="tag warn">differs from artifact</span>`}</td></tr>`).join("")}
      </tbody></table>`;
    const diff = [];
    for (const cand of cands) {
      const now = assigned.get(cand.variant_key);
      if (now.tier !== cand.tier) diff.push({ key: cand.variant_key, old: cand.tier, now: now.tier, r: now });
    }
    $("#pl_out_diff").innerHTML = diff.length
      ? `<table><thead><tr><th>variant</th><th>artifact tier</th><th>→ under your rule</th><th>why</th></tr></thead><tbody>
        ${diff.slice(0, 12).map((x) => `<tr><td class="small mono">${esc(x.key)}</td><td>${tierBadge(x.old)}</td><td>${tierBadge(x.now)}</td>
          <td class="wrap small">${esc(gateExplain(x.r))}</td></tr>`).join("")}</tbody></table>
          ${diff.length > 12 ? `<p class="small dim">…and ${nfmt(diff.length - 12)} more.</p>` : ""}`
      : `<p class="small dim">no candidate changes tier under this rule.</p>`;
  };
  for (const [id, out] of [["#pl_nm", "#pl_nm_v"], ["#pl_tpm", "#pl_tpm_v"], ["#pl_vaf", "#pl_vaf_v"]]) {
    $(id).addEventListener("input", () => { $(out).textContent = $(id).value; applyPlay(); });
  }
  applyPlay();
}

/* ---------------- L4: TESLA benchmark ---------------- */

function viewTesla() {
  const t = state.tesla;
  if (!t) {
    app.innerHTML = `<h1>Ranking-rule audit — TESLA ${layerTag("l4", "L4 experimental")}</h1>
      <div class="notice"><span class="ttl">Pending.</span> No exported benchmark artifact was found. Until the
      evaluation harness has been run, TESLA validation stays marked pending — this interface never displays a
      number it cannot trace to a hashed artifact.</div>`;
    return;
  }
  const m = t.metrics;
  const rows = rowsOf(t.validated_predictions);
  app.innerHTML = `
    <h1>Ranking-rule audit on the TESLA benchmark ${layerTag("l4", "L4 experimental")}</h1>
    <p class="sub">TESLA is a community benchmark in which predicted peptide–MHC candidates were
    <em>experimentally</em> measured. NEOFORGE uses it for exactly one question: does the priority ranking rule
    (order by predicted affinity) place experimentally-validated pMHC near the top more often than chance?
    That is cohort-level evidence about a rule. It does not validate the HCC1395 case, and it does not predict
    immunogenicity of anything shown in this demo.</p>
    <p class="small dim mono">${esc(t.source)}</p>

    <div class="card">
      <h3 style="margin-top:0">Metrics</h3>
      <div class="stats">
        <div class="stat"><div class="v">${num(m.ttif_top20pct, 3)}</div><div class="k">TTIF top 20% · random ${num(m.ttif_random_expectation, 2)} · oracle ${num(m.ttif_oracle, 2)}</div></div>
        <div class="stat"><div class="v">${num(m.auprc, 3)}</div><div class="k">AUPRC · prevalence ${num(m.auprc_random_expectation, 3)} · oracle ${num(m.auprc_oracle, 2)}</div></div>
        <div class="stat"><div class="v">${nfmt(m.median_rank_of_validated)}</div><div class="k">median rank of validated · of ${nfmt(m.n_candidates)} pMHC</div></div>
        <div class="stat"><div class="v">${nfmt(m.validated_in_top20)} / ${nfmt(m.n_validated)}</div><div class="k">validated in top 20% · random expect ${num(m.validated_in_top20_random_expectation, 1)}</div></div>
      </div>
    </div>

    <div class="${t.rank_rule_falsified ? "notice red" : "notice good"}">
      <span class="ttl">Falsification verdict: ${t.rank_rule_falsified ? "the ranking rule IS falsified on this benchmark" : "the ranking rule is NOT falsified on this benchmark"}.</span><br>
      <span class="mono small">pre-registered clause: ${esc(t.parameters.falsification_clause)}</span><br>
      <span class="small">Had the measured ${num(m.ttif_top20pct, 3)} been ≤ ${num(m.ttif_random_expectation, 1)}, the priority design would have been
      declared falsified and this page would say so.</span>
    </div>

    <div class="card">
      <h3 style="margin-top:0">Methodology — exactly what was run</h3>
      <dl class="kv">
        <dt>predictor</dt><dd>${esc(t.parameters.predictor)} (MHCflurry ${esc(t.parameters.mhcflurry_version)})</dd>
        <dt>labels file</dt><dd class="wrap">${esc(t.parameters.labels)}</dd>
        <dt>labels layer</dt><dd>${layerTag("l4", "L4 measured binding")} — experimental, third-party</dd>
        <dt>allele normalisation</dt><dd>${esc(t.parameters.allele_normalisation)}</dd>
        <dt>ranking</dt><dd>${esc(t.parameters.ranking)}</dd>
        <dt>scale</dt><dd>${nfmt(m.n_candidates)} pMHC · ${m.n_patients} patients · ${m.n_alleles} alleles</dd>
      </dl>
    </div>

    <div class="two-col">
      <div class="card">
        <h3 style="margin-top:0">Fraction of validated pMHC found in the top k% of the ranking</h3>
        <p class="small dim">Solid line: the published affinity ranking. Dashed: analytic random expectation (k/100).
        Every point above the dashed line is enrichment the rule earns; the k=20% point is the TTIF headline.</p>
        ${curveSVG(m.rank_at_k_curve)}
      </div>
      <div class="card">
        <h3 style="margin-top:0">Rule ablation — same labels, different ordering criteria</h3>
        ${t.rule_ablation && t.rule_ablation.entries.length ? `
        <table><thead><tr><th>ordering criterion</th><th class="num">TTIF top 20%</th>
          <th class="num">AUPRC</th><th class="num">median rank</th></tr></thead><tbody>
          ${t.rule_ablation.entries.map((e) => `
            <tr><td class="wrap small">${esc(e.criterion)}${e.criterion.includes("published") ? " ← rule under test" : ""}</td>
            <td class="num">${num(e.ttif_top20pct, 3)}</td>
            <td class="num">${num(e.auprc, 3)}</td>
            <td class="num">${num(e.median_rank_of_validated)}</td></tr>`).join("")}
        </tbody></table>
        <p class="small dim">${esc(t.rule_ablation.note)}</p>` :
        `<div class="notice small"><span class="ttl">Ablation unavailable:</span>
         ${esc(t.rule_ablation ? t.rule_ablation.note : "not present in this export")} — recorded, not imputed.</div>`}
      </div>
    </div>

    <div class="card">
      <h3 style="margin-top:0">The ${nfmt(rows.length)} experimentally-validated pMHC</h3>
      <p class="small dim">Predicted affinity in this table is still a prediction ${layerTag("l3")};
      only <span class="mono">measured nM</span> is experimental ${layerTag("l4")}. The two columns are never merged.</p>
      <div style="max-height:420px;overflow-y:auto;border:1px solid var(--line-soft);border-radius:8px">
      <table><thead><tr><th>peptide</th><th>allele</th><th class="num">rank</th>
        <th class="num">pred. affinity nM (L3)</th><th class="num">measured nM (L4)</th></tr></thead><tbody>
        ${rows.map((r) => `<tr><td class="mono">${esc(r.peptide)}</td><td class="small">${esc(r.allele)}</td>
          <td class="num">${nfmt(r.rank)}</td>
          <td class="num">${num(r.mhcflurry_affinity_nm, 1)}</td>
          <td class="num">${num(r.measured_binding_affinity_nm, 1)}</td></tr>`).join("")}
      </tbody></table></div>
    </div>

    <div class="card">
      <h3 style="margin-top:0">Non-claims for this benchmark</h3>
      <ul class="small dim">${t.non_claims.map((c) => `<li>${esc(c)}</li>`).join("")}</ul>
    </div>`;
}

/* ---------------- workbench ---------------- */

function viewWorkbench() {
  const d = state.data;
  const b = d.case_bundle;
  const stages = Object.entries(d.pipeline_stage_records);
  app.innerHTML = `
    ${caseHeader(d, "Analysis workbench", "The same data the explorer shows — behind it: files, checksums, parameters, versions, artifact hashes, timings, uncertainty.")}

    <details class="full" open><summary>Case bundle — inputs pinned by SHA-256</summary><div>
      <p class="small dim">${esc(b.description)}</p>
      <dl class="kv">
        <dt>genome build</dt><dd>${esc(b.genome_build || "—")}</dd>
        ${b.pairing ? `<dt>pairing</dt><dd class="wrap">${esc(typeof b.pairing === "string" ? b.pairing : JSON.stringify(b.pairing))}</dd>` : ""}
      </dl>
      <table><thead><tr><th>input</th><th>path / role</th><th>sha256</th></tr></thead><tbody>
        ${(d.ingest.inputs || []).map((i) => `<tr><td>${esc(i.name)}</td>
          <td class="wrap"><span class="mono small">${esc(i.path || "")}</span><br><span class="dim small">${esc(i.role || "")}</span></td>
          <td class="sha">${esc(i.sha256 || "")}</td></tr>`).join("")}
      </tbody></table>
      ${b.expression ? `<dl class="kv"><dt>expression</dt><dd class="wrap">${esc(b.expression.matrix || (b.expression.patient_sample ? "same-patient GDC STAR counts: " + b.expression.patient_sample : JSON.stringify(b.expression)))}</dd></dl>` : ""}
    </div></details>

    <details class="full"><summary>Stage artifacts — the content-addressed chain</summary><div>
      <table><thead><tr><th>stage</th><th>artifact id</th><th>sha256</th><th class="num">wall (s)</th></tr></thead><tbody>
        ${stages.map(([s, r]) => `<tr><td>${esc(s)}</td><td class="mono small">${esc(r.artifact_id)}</td>
          <td class="sha">${esc(r.sha256 || "")}</td><td class="num">${r.wall_seconds ?? "—"}</td></tr>`).join("")}
      </tbody></table>
      <p class="small dim">Re-running the pipeline reproduces these artifact ids exactly (determinism is covered by the
      test suite). The exported projections on every screen come from these same artifacts — no side-channel
      presentation data exists.</p>
    </div></details>

    <details class="full"><summary>Stage parameters</summary><div>
      ${["variants", "translate", "present", "prioritize", "report"].map((s) => {
        const blk = d[s];
        const rec = d.pipeline_stage_records[s];
        const params = (blk && blk.parameters) || (rec && rec.parameters) || null;
        return params ? `<h3>${esc(s)}</h3><dl class="kv">
          ${Object.entries(params).map(([k, v]) => `<dt>${esc(k)}</dt><dd class="wrap">${esc(typeof v === "object" ? JSON.stringify(v) : String(v))}</dd>`).join("")}
        </dl>` : "";
      }).join("")}
    </div></details>

    <details class="full"><summary>Validation checks</summary><div>
      <table><tbody>${(d.ingest.checks || []).map((c) => `<tr><td style="width:56px">${statusTag(c.status)}</td>
        <td class="mono" style="width:340px">${esc(c.id)}</td><td class="wrap">${esc(c.detail)}</td></tr>`).join("")}</tbody></table>
    </div></details>

    ${d.report ? `
    <details class="full"><summary>Uncertainty &amp; limitation register</summary><div>
      <table><thead><tr><th>status</th><th>source</th><th>detail</th></tr></thead><tbody>
        ${d.report.uncertainty_register.map((u) => `<tr><td>${statusTag(u.status)}</td>
          <td class="mono small">${esc(u.source)}</td><td class="wrap small">${esc(u.detail)}</td></tr>`).join("")}
      </tbody></table>
    </div></details>

    <details class="full"><summary>What science would have to confirm — and the non-claims</summary><div>
      <ul class="small">${(d.report.what_science_would_confirm || []).map((c) => `<li>${esc(c)}</li>`).join("")}</ul>
      <h3>Non-claims</h3>
      <ul class="small dim">${(d.report.non_claims || []).map((c) => `<li>${esc(c)}</li>`).join("")}</ul>
    </div></details>

    <details class="full"><summary>Evidence ledger — which layer supplied what</summary><div>
      ${["L1_measured", "L2_derived", "L3_predicted", "L4_experimental"].map((k) => {
        const led = d.report.evidence_ledger[k];
        return led ? `<h3>${layerTag(k.slice(0, 2).toLowerCase(), k.replace("_", " ").toLowerCase())}</h3>
          <dl class="kv">${Object.entries(led).filter(([kk]) => kk !== "kind")
            .map(([kk, v]) => `<dt>${esc(kk)}</dt><dd class="wrap">${esc(typeof v === "object" ? JSON.stringify(v) : String(v))}</dd>`).join("")}</dl>` : "";
      }).join("")}
      ${d.report.evidence_ledger.tesla_cohort_reference ? `<div class="notice good small">
        <span class="ttl">L4 cohort reference embedded:</span>
        ${esc(JSON.stringify(d.report.evidence_ledger.tesla_cohort_reference))}</div>`
        : `<div class="notice small"><span class="ttl">L4 experimental:</span> no benchmark artifact present in this
        export — the ledger records exactly that instead of omitting the layer silently.</div>`}
    </div></details>` : ""}

    ${state.tesla ? `<details class="full"><summary>Benchmark artifact (L4 cohort reference)</summary><div>
      <p class="small dim mono">${esc(state.tesla.source)}</p>
      <p class="small dim">Referenced by the report stage as cohort-level ranking-rule evidence only — never as
      validation of this case. See <a href="#/tesla">the benchmark page</a>.</p>
    </div></details>` : ""}

    <details class="full"><summary>Export provenance — the bridge module</summary><div>
      <p class="small dim mono">${esc(d.exported_by)}</p>
      <p class="small dim">This interface renders only these exported projections. The projection module
      (<span class="mono">neoforge.export_ui</span>) contains no scientific computation, is deliberately not a
      pipeline stage, and adds no claim and no number.</p>
    </div></details>`;
}

/* ---------------- research report ---------------- */

function renderMd(src) {
  const inline = (s) => s
    .replace(/\*\*([^*]+)\*\*/g, "<strong>$1</strong>")
    .replace(/`([^`]+)`/g, "<code>$1</code>")
    .replace(/\*([^*]+)\*/g, "<em>$1</em>");
  const lines = esc(src).split("\n");
  let html = "", listOpen = false, inTable = false;
  const flushList = () => { if (listOpen) { html += "</ul>"; listOpen = false; } };
  const closeTable = () => { if (inTable) { html += "</tbody></table>"; inTable = false; } };
  for (const raw of lines) {
    const line = raw.replace(/\s+$/, "");
    if (/^\s*$/.test(line)) { flushList(); closeTable(); continue; }
    if (/^#{1,4} /.test(line)) {
      flushList(); closeTable();
      const lvl = Math.min(line.match(/^#+/)[0].length, 4);
      html += `<h${lvl}>${inline(line.replace(/^#+ /, ""))}</h${lvl}>`;
      continue;
    }
    if (/^\s*---+\s*$/.test(line)) { flushList(); closeTable(); html += "<hr>"; continue; }
    if (/^\s*-\s+/.test(line)) {
      closeTable();
      if (!listOpen) { html += "<ul>"; listOpen = true; }
      html += `<li>${inline(line.replace(/^\s*-\s+/, ""))}</li>`;
      continue;
    }
    if (line.includes("|")) {
      const cells = line.split("|").map((c) => c.trim());
      const isSeparator = cells.every((c) => c === "" || /^:?-{2,}:?$/.test(c));
      if (isSeparator) continue;
      if (!inTable) {
        flushList();
        html = html.replace(/<h([1-4])>([^<]*)<\/h\1>\s*$/, '<h$1>$2</h$1>'); // keep heading text intact
        html += "<table><tbody>";
        inTable = true;
      }
      const trimmed = cells.slice(cells[0] === "" ? 1 : 0, cells[cells.length - 1] === "" ? cells.length - 1 : cells.length);
      html += `<tr>${trimmed.map((c) => `<td>${inline(c)}</td>`).join("")}</tr>`;
      continue;
    }
    flushList(); closeTable();
    html += `<p>${inline(line)}</p>`;
  }
  flushList(); closeTable();
  return `<div class="md">${html}</div>`;
}

function viewReport() {
  const d = state.data;
  if (!d.report) {
    app.innerHTML = caseHeader(d, "Research report", "") +
      `<div class="notice"><span class="ttl">No research report for this case.</span> The report stage reads the
      prioritisation output, which requires a same-patient HLA genotype. This case's chain honestly stops at L2
      (peptide derivation). The workbench still exposes its full provenance, and
      <a href="#/mutations">the mutation explorer</a> remains fully available.</div>`;
    return;
  }
  const r = d.report;
  app.innerHTML = `
    ${caseHeader(d, "Research report", "An exportable, evidence-layered summary for qualified researchers — explicitly not a treatment recommendation.")}
    <div class="controls">
      <button id="rep_dl">download .md</button>
      <button id="rep_copy">copy markdown</button>
    </div>
    <div class="card">${renderMd(r.research_report_md)}</div>
    <div class="card">
      <h3 style="margin-top:0">Evidence ledger — which layer supplied what</h3>
      ${["L1_measured", "L2_derived", "L3_predicted", "L4_experimental"].map((k) => {
        const led = r.evidence_ledger[k];
        return led ? `<details class="full"><summary>${layerTag(k.slice(0, 2).toLowerCase(), k.replace("_", " "))}</summary><div>
          <dl class="kv">${Object.entries(led).filter(([kk]) => kk !== "kind")
            .map(([kk, v]) => `<dt>${esc(kk)}</dt><dd class="wrap">${esc(typeof v === "object" ? JSON.stringify(v) : String(v))}</dd>`).join("")}</dl>
        </div></details>` : "";
      }).join("")}
      ${r.evidence_ledger.tesla_cohort_reference ? `
        <div class="notice good small"><span class="ttl">L4 cohort reference embedded.</span>
        Ranking rule audited on TESLA: TTIF top-20% ${num(r.evidence_ledger.tesla_cohort_reference.ttif_top20pct, 3)}
        (random expectation ${num(r.evidence_ledger.tesla_cohort_reference.ttif_random_expectation, 2)}),
        AUPRC ${num(r.evidence_ledger.tesla_cohort_reference.auprc, 3)},
        rule ${r.evidence_ledger.tesla_cohort_reference.rank_rule_falsified ? "FALSIFIED" : "not falsified"}.
        ${esc(r.evidence_ledger.tesla_cohort_reference.scope_note || "")}</div>` :
        `<div class="notice small"><span class="ttl">L4 experimental:</span> no benchmark artifact present in this
        export — the ledger records exactly that rather than dropping the layer silently.</div>`}
    </div>`;

  $("#rep_dl").addEventListener("click", () => {
    const url = URL.createObjectURL(new Blob([r.research_report_md], { type: "text/markdown" }));
    const a = document.createElement("a");
    a.href = url;
    a.download = `neoforge_research_report_${d.case_bundle.case_id}.md`;
    a.click();
    URL.revokeObjectURL(url);
  });
  $("#rep_copy").addEventListener("click", async (e) => {
    try { await navigator.clipboard.writeText(r.research_report_md); e.target.textContent = "copied ✓"; }
    catch { e.target.textContent = "clipboard blocked"; }
    setTimeout(() => { e.target.textContent = "copy markdown"; }, 1600);
  });
}

/* ---------------- construct (design layer) ---------------- */

function constructMapSVG(c) {
  const segments = c.nucleotide_map;
  const aaTotal = segments.reduce((s, seg) => s + (seg.aa ? seg.aa.length : 0), 0);
  if (!aaTotal) return "";
  const W = 940, H = 92, pad = 10, barY = 34, barH = 26;
  const scale = (W - 2 * pad) / aaTotal;
  const flagByJunction = {};
  for (const r of c.junction_audit.final_rows) flagByJunction[r.junction] = flagByJunction[r.junction] || (r.flagged ? 1 : 0);
  const jColors = { placeholder: ["none", "var(--faint)", "6 3"],
                    design_assumption: ["var(--accent-soft)", "var(--accent)", ""],
                    linker: ["#58a6ff33", "var(--l2)", ""],
                    selected_epitope: ["var(--accent-soft)", "var(--accent)", ""] };
  const epitopeEnds = [];
  let aa = 0;
  let svg = "", chips = [];
  for (const seg of segments) {
    const isAA = !!seg.aa;
    const w = isAA ? Math.max(3, seg.aa.length * scale) : 14;
    const [fill, stroke, dash] = jColors[seg.kind] || ["none", "var(--line)", ""];
    const x = isAA ? pad + aa * scale : pad + Math.min(aa * scale, W - pad - w);
    if (seg.kind === "selected_epitope") epitopeEnds.push({ end: aa + seg.aa.length, flag: null });
    if (isAA) aa += seg.aa.length;
    svg += `<rect x="${x.toFixed(1)}" y="${barY}" width="${w.toFixed(1)}" height="${barH}" fill="${fill}" stroke="${stroke}" stroke-width="1" ${dash ? `stroke-dasharray="${dash}"` : ""} rx="2"><title>${esc(seg.name)}${seg.aa ? ` \u00b7 ${seg.aa.length} aa` : ""} \u00b7 ${esc(seg.kind)}${seg.note ? " \u00b7 " + esc(seg.note) : ""}</title></rect>`;
    if (seg.kind === "selected_epitope") {
      const idx = epitopeEnds.length;
      if (w > 26) svg += `<text class="axis" x="${(x + w / 2).toFixed(1)}" y="${barY + 17}" text-anchor="middle">E${idx}</text>`;
      chips.push(`<span class="tag ok" title="${esc(seg.allele)} \u00b7 ${num(seg.mutant_affinity_nm, 1)} nM predicted">E${idx} ${esc(seg.aa)}</span>`);
    } else if (w > 18) {
      svg += `<text class="axis" x="${(x + w / 2).toFixed(1)}" y="${barY + 44}" text-anchor="middle">${esc(seg.name.split(":")[0].split(" ")[0])}</text>`;
    }
  }
  for (let j = 0; j < epitopeEnds.length - 1; j++) {
    const bx = pad + (epitopeEnds[j].end + 1.5) * scale;
    const flagged = (flagByJunction[j] || 0) > 0;
    svg += flagged
      ? `<polygon points="${bx.toFixed(1)},${barY - 8} ${(+bx - 4).toFixed(1)},${barY - 1} ${(+bx + 4).toFixed(1)},${barY - 1}" fill="var(--bad)"><title>junction ${j}: creates newly flagged strong binders across at least one allele</title></polygon>`
      : `<line x1="${bx.toFixed(1)}" y1="${barY - 6}" x2="${bx.toFixed(1)}" y2="${barY - 1}" stroke="var(--faint)"/><title>junction ${j}: no newly flagged strong binder</title>`;
  }
  return `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="schematic construct map"
    style="max-width:${W}px;width:100%">
    <line x1="${pad}" y1="${barY + barH}" x2="${W - pad}" y2="${barY + barH}" stroke="var(--line)"/>
    <text class="axis" x="${pad}" y="12">schematic research-design map \u2014 block lengths proportional to residue count; NOT a predicted molecular structure</text>
    ${svg}</svg>
    <p class="small" style="line-height:2.1">
      <span class="tag l4" title="declared design assumption, not a validation">design assumption</span>
      <span class="tag l2" title="MHC-I AAY linkers between epitopes">linker ${esc(c.nucleotide_stats.linker)}</span>
      placeholders (dashed): 5'UTR, 3'UTR, poly(A) \u2014 labelled, not simulated.
      &nbsp;${chips.join(" ")}</p>`;
}

function structureStripSVG(s) {
  if (!s || !s.propensity_table) return "";
  const labels = s.propensity_table.label;
  const aa = s.propensity_table.aa;
  const n = labels.length;
  if (!n) return "";
  const W = 940, pad = 10, H = 58, barY = 22, barH = 18;
  const scale = (W - 2 * pad) / n;
  const colors = { H: "var(--accent)", E: "#58a6ff", C: "var(--faint)" };
  const names = { H: "alpha-helix propensity", E: "beta-sheet propensity", C: "coil (reference 1.0)" };
  let rects = "";
  for (let i = 0; i < n; i++) {
    rects += `<rect x="${(pad + i * scale).toFixed(2)}" y="${barY}" width="${Math.max(scale, 0.8).toFixed(2)}" height="${barH}" fill="${colors[labels[i]]}"><title>residue ${i + 1} ${esc(aa[i])} \u00b7 ${labels[i]} (${esc(names[labels[i]] || labels[i])})</title></rect>`;
  }
  const counts = s.summary.label_counts;
  const runN = s.summary.helix_runs_ge3.length + s.summary.sheet_runs_ge3.length;
  return `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="secondary-structure propensity strip"
    style="max-width:${W}px;width:100%">
    <text class="axis" x="${pad}" y="12">residue 1..${nfmt(n)} \u00b7 one bar per residue \u00b7 NOT a 3-D model</text>
    ${rects}</svg>
    <p class="small">
      <span class="tag ok">H ${nfmt(counts["H"] || 0)}</span>
      <span class="tag l2">E ${nfmt(counts["E"] || 0)}</span>
      <span class="tag missing">C ${nfmt(counts["C"] || 0)}</span>
      <span class="tag warn" title="runs of \u22653 consecutive residues with the same label">${nfmt(runN)} runs \u22653</span>
      &nbsp;windowed mean of pinned Chou\u2013Fasman (1978) parameters, argmax over {H, E, C};
      ${esc(s.parameters.chou_fasman.file)} \u00b7 sha256 ${shortSha(s.parameters.chou_fasman.sha256)}</p>`;
}

function viewConstruct() {
  const d = state.data;
  const c = d.construct;
  if (!c) {
    app.innerHTML = caseHeader(d, "Construct (design)", "Why is there no construct design for this case?") +
      `<div class="notice red"><span class="ttl">No construct artifact in this case's export.</span><br>
      The construct stage is mechanically gated: it requires the full L1\u2013L3 chain (including HLA-resolved
      presentation predictions) and at least one Tier A candidate that survives every documented selection rule.
      This case's chain stops before those gates \u2014 and a refused stage must not look like missing polish;
      it <em>is</em> the honesty mechanism. ${stageSourceLine(d, "prioritize")}</div>`;
    return;
  }
  const cs = c.case_summary, st = c.nucleotide_stats;
  const flags = c.junction_audit.final_flagged_junction_peptides;
  const byJunction = {};
  for (const r of c.junction_audit.final_rows) {
    byJunction[r.junction] = byJunction[r.junction] || { rows: 0, flagged: 0, peptides: new Set() };
    const b = byJunction[r.junction];
    b.rows++; if (r.flagged) b.flagged++;
    b.peptides.add(r.peptide);
  }
  const jEntries = Object.keys(byJunction).map(Number).sort((a, b) => a - b)
    .map((j) => ({ label: `J${j + 1}`, value: byJunction[j].flagged,
                   color: byJunction[j].flagged ? "var(--bad)" : "var(--faint)" }));
  const flaggedRows = c.junction_audit.final_rows.filter((r) => r.flagged);
  const nCheck = (k, good) => `<span class="tag ${good ? "ok" : "bad"}">${esc(k)}</span>`;
  app.innerHTML = `
    ${caseHeader(d, "Construct (design)", "From the prioritised tier-A shortlist to an in-silico polyepitope research design \u2014 what would a vaccine-research construct look like, reductively, and why every choice is ledgered.")}
    <div class="notice red"><span class="ttl">${esc(c.stamp)}</span><br>
      ${esc(c.non_claims.join(" "))}</div>
    ${stageSourceLine(d, "construct")}

    <div class="card">
      <h3 style="margin-top:0">Construct map \u2014 ${nfmt(st.protein_length_aa)} aa encoded by ${nfmt(st.cds_length_nt)} nt</h3>
      ${constructMapSVG(c)}
    </div>

    ${d.structure ? `<div class="card">
      <h3 style="margin-top:0">Secondary-structure propensity \u2014 sequence-derived, not a 3-D model</h3>
      <p class="small dim">${esc(d.structure.stamp)} ${esc(d.structure.non_claims.join(" "))}
      <br>validation: <span class="tag missing">${esc(d.structure.validation_status)}</span>
      ${esc(d.structure.validation_note)}</p>
      ${structureStripSVG(d.structure)}
      ${stageSourceLine(d, "structure")}
    </div>` : ""}


    <div class="card">
      <h3 style="margin-top:0">Selection ledger \u2014 ${nfmt(cs.epitopes_selected)} included \u00b7 ${nfmt(cs.epitopes_excluded_ledgered)} ledgered exclusions</h3>
      <div class="stats">
        <div class="stat"><div class="v">${nfmt(cs.epitopes_selected)}</div><div class="k">epitopes included \u00b7 tier-A pool ${nfmt(cs.tier_a_candidates)} \u00b7 cap ${nfmt(c.parameters.max_epitopes)}</div></div>
        <div class="stat"><div class="v">${nfmt(cs.alleles_covered.length)} / ${nfmt(cs.case_alleles.length)}</div><div class="k">case HLA alleles covered by the selection</div></div>
        <div class="stat"><div class="v">${nfmt(cs.genes_covered.length)}</div><div class="k">genes covered (diversity preference)</div></div>
        <div class="stat"><div class="v">${nfmt(c.exclusion_reason_counts["no_mutant_vs_wild_type_discrimination"] || 0)}</div><div class="k">dropped: mutant\u2260wild-type discrimination absent</div></div>
      </div>
      <table><thead><tr><th>variant</th><th>peptide</th><th>allele</th><th>gene</th>
        <th class="num">mutant nM (L3)</th><th class="num">ref nM (L3)</th><th>selection reason</th></tr></thead><tbody>
        ${c.selection_ledger_included.map((e) => `
          <tr><td class="mono"><a href="#/protein?variant=${encodeURIComponent(e.variant_key)}">${esc(e.variant_key)}</a></td>
          <td class="mono">${esc(e.peptide)}</td><td class="small">${esc(e.allele)}</td>
          <td class="small mono">${esc(e.gene_id)}</td>
          <td class="num">${num(e.mutant_affinity_nm, 1)}</td>
          <td class="num">${num(e.reference_affinity_nm, 1)}</td>
          <td class="small">${esc(e.selection_reason)}</td></tr>`).join("")}
      </tbody></table>
      <h3>Exclusions \u2014 every mutation leaves the pool through a documented door</h3>
      <p class="small">${Object.entries(c.exclusion_reason_counts).map(([r, n]) =>
        `<span class="tag ${r === "not_tier_A" ? "missing" : "warn"}" title="${esc(r)}">${nfmt(n)} \u00d7 ${esc(r)}</span>`).join(" ")}</p>
      <details class="full"><summary>inspect all ${nfmt(c.selection_ledger_excluded.length)} ledgered exclusions</summary><div style="max-height:320px;overflow-y:auto">
      <table><thead><tr><th>variant</th><th>reason</th><th>detail</th></tr></thead><tbody>
        ${c.selection_ledger_excluded.map((e) => `
          <tr><td class="mono small">${esc(e.variant_key)}</td><td class="small">${esc(e.reason)}</td>
          <td class="small dim">${esc(e.detail || "")}</td></tr>`).join("")}
      </tbody></table></div></details>
    </div>

    <div class="two-col">
      <div class="card">
        <h3 style="margin-top:0">Junction audit \u2014 what the joins might create</h3>
        <p class="small dim">Every 8\u201311-mer spanning an epitope\u2013${esc(c.parameters.linker)} join was scored
        ${layerTag("l3", "L3 predicted")} with ${esc(c.junction_audit.predictor.software)} against every case
        allele (<span class="mono">${nfmt(c.junction_audit.final_rows.length)}</span> rows). A new strong binder
        (${num(c.parameters.junction_flag_threshold_nm, 0)} nM) created by a join is flagged, never averaged away.</p>
        <div class="stats">
          <div class="stat"><div class="v" style="color:${flags ? "var(--bad)" : "var(--l1)"}">${nfmt(flags)}</div><div class="k">final flagged junction peptides \u00b7 initial ordering had ${nfmt(c.junction_audit.history[0].flagged_junction_peptides)}</div></div>
        </div>
        ${jEntries.length ? `<p class="small dim">flagged peptides per junction:</p>${barChartSVG(jEntries, { width: 470, labelW: 34 })}` : ""}
        <details class="full"><summary>audit history (${c.junction_audit.history.length} steps, accepted orderings must be monotone)</summary><div>
          <table><thead><tr><th>action</th><th class="num">flagged</th><th>accepted</th></tr></thead><tbody>
          ${c.junction_audit.history.map((h) => `<tr><td class="small">${esc(h.action)}</td>
            <td class="num">${nfmt(h.flagged_junction_peptides)}</td>
            <td>${h.accepted ? '<span class="tag ok">yes</span>' : '<span class="tag missing">no</span>'}</td></tr>`).join("")}
          </tbody></table></div></details>
        <details class="full"><summary>all flagged rows (${nfmt(flaggedRows.length)})</summary><div style="max-height:280px;overflow-y:auto;border:1px solid var(--line-soft);border-radius:8px">
          <table><thead><tr><th class="num">junction</th><th>peptide</th><th class="num">len</th><th>allele</th>
          <th class="num">pred. nM</th></tr></thead><tbody>
          ${flaggedRows.map((r) => `<tr><td class="num">J${r.junction + 1}</td><td class="mono">${esc(r.peptide)}</td>
            <td class="num">${r.length}</td><td class="small">${esc(r.allele)}</td>
            <td class="num">${num(r.mhcflurry_affinity_nm, 1)}</td></tr>`).join("")}
          </tbody></table></div></details>
      </div>
      <div class="card">
        <h3 style="margin-top:0">Back-transcription \u2014 codons, checked, not asserted</h3>
        <p class="small">${nCheck("translation matches design", c.sequence_checks.translation_matches_design)}
          ${nCheck("0 internal stop codons", c.sequence_checks.internal_stop_codons === 0)}
          ${nCheck("starts with ATG", c.sequence_checks.cds_starts_with_atg)}
          <span class="tag missing" title="UTR / poly(A) are labelled placeholders">UTRs + poly(A) placeholders</span>
          <span class="tag missing" title="N1-methylpseudouridine, modified nucleosides, formulation \u2014 not simulated">mRNA chemistry not simulated</span></p>
        <dl class="kv">
          <dt>CDS length</dt><dd>${nfmt(st.cds_length_nt)} nt (${nfmt(st.protein_length_aa)} aa + ${esc(c.sequence_checks.stop_codon)} stop)</dd>
          <dt>GC content</dt><dd>${num(st.gc_percent, 2)}%</dd>
          <dt>codon table</dt><dd class="wrap small">${esc(c.parameters.codon_table.derivation)} \u00b7 sha256 <span class="mono">${shortSha(c.parameters.codon_table.sha256)}</span></dd>
          <dt>homopolymer mitigation</dt><dd>${esc(c.reverse_translation.codon_choices_summary)} (threshold ${st.codon_run_threshold ?? c.parameters.max_homopolymer_run})</dd>
          <dt>homopolymer runs \u2265 5</dt><dd>${st.homopolymer_runs_ge5.length ? st.homopolymer_runs_ge5.map((r) => `${esc(r.base)}\u00d7${r.length} @${r.start}`).join(", ") : "none"}</dd>
          <dt>restriction sites</dt><dd>${Object.keys(st.restriction_sites).length
            ? Object.entries(st.restriction_sites).map(([k, v]) => `${esc(k)} ${esc(v.site)} \u00d7${v.positions_0based.length}`).join(", ")
            : "none in the ${nfmt(12)}-enzyme documented panel"}</dd>
          <dt>linkers</dt><dd>${nfmt(st.linker_count)} \u00d7 ${esc(st.linker)} (published MHC-I convention)</dd>
          <dt>signal peptide</dt><dd class="wrap small">${esc(c.parameters.signal_peptide)} \u2014 ${esc(c.parameters.signal_peptide_source)}</dd>
          <dt>trafficking tag</dt><dd class="wrap small">${esc(c.parameters.trafficking_tag)} \u2014 ${esc(c.parameters.trafficking_tag_source)}</dd>
          <dt>FASTA export</dt><dd><a href="data/${esc(c.fasta_file)}" download="${esc(c.fasta_file)}"><button class="controls">download .fasta</button></a>
            <span class="mono small" title="sha256 of the served copy (byte-identical to the artifact-store file)">${esc(shortSha(c.fasta_sha256))}</span></dd>
        </dl>
      </div>
    </div>

    <div class="two-col">
      <div class="card">
        <h3 style="margin-top:0">Design assumptions (declared, not validated)</h3>
        <ul class="small dim">${c.design_assumptions.map((a) => `<li>${esc(a)}</li>`).join("")}</ul>
      </div>
      <div class="card">
        <h3 style="margin-top:0">Uncertainty register</h3>
        <ul class="small dim">${c.uncertainty_register.map((a) => `<li>${esc(a)}</li>`).join("")}</ul>
        <p class="small dim">Inputs, hashed: ${Object.entries(c.inputs).map(([k, v]) =>
          `${esc(k)} <span class="mono">${esc(v.artifact_id)}</span>`).join(" \u00b7 ")}</p>
      </div>
    </div>`;
}

/* ---------------- rosie pathway mapping ---------------- */

const ROSIE_STEPS = [
  { step: "1 \u00b7 Tumour + normal-tissue sequencing", target: "Stage 1 \u2014 ingest (L1 measured)",
    does: "Consumes public, licensed sequencing files; validates tumour\u2013normal pairing, depth, TVAF coverage and HLA typing, with per-input SHA-256 provenance.",
    not: "NEOFORGE does not sequence anything \u2014 no sample acquisition, no novel patient data." },
  { step: "2 \u00b7 Mutation analysis (tumour vs normal)", target: "Stage 2 \u2014 variants (L1)",
    does: "Tumour-vs-normal SNV comparison with a conservation ledger; measured TVAF and depth carried through; patient-barcode gates for same-patient integrity.",
    not: "No de novo variant discovery beyond documented filters; INDELs are recorded out of scope (v0.1), not hidden." },
  { step: "3 \u00b7 Neoantigen investigation (derive fragments)", target: "Stage 3 \u2014 translate (L2 derived)",
    does: "Derives mutant vs reference peptide windows from genome + GENCODE; every derived window verified against official protein translations.",
    not: "Not a structure model; no proteasomal-processing or transport prediction beyond what MHCflurry provides." },
  { step: "4 \u00b7 HLA-restricted presentation assessment", target: "Stage 4 \u2014 present (L3 predicted)",
    does: "Per-allele MHCflurry binding + presentation predictions, mutant vs reference windows, wrapper-only use with versioned provenance.",
    not: "Predicted binding is NOT predicted immunogenicity; unsupported alleles are refused, never imputed." },
  { step: "5 \u00b7 Candidate assessment & prioritisation", target: "Stage 5 \u2014 prioritise + L4 TESLA audit",
    does: "Explainable tiering on published thresholds (500 nM, TPM, VAF) \u2014 no composite score; ranking rule audited against TESLA experimental labels (cohort-level).",
    not: "No candidate-level immunogenicity label exists for this case's variants; the TESLA audit validates a rule, not a patient." },
  { step: "6 \u00b7 Protein modelling", target: "\u2014 not implemented in this pass",
    does: "Declared deferred: conceptual molecular illustrations are labelled as conceptual in the interface.",
    not: "No structure prediction is presented as performed; no mislabelled protein models." },
  { step: "7 \u00b7 Construct design (mRNA)", target: "Stage 7 \u2014 construct (design layer)",
    does: "In-silico polyepitope research design from tier-A candidates: ledgered selection, MHC-I linkers, junction audit, pinned-codon back-translation, hashed artifact.",
    not: "A research design stamp \u2014 not manufacturable, not validated, not a vaccine; mRNA chemistry not simulated." },
  { step: "8 \u00b7 Scientific collaboration & review", target: "Provided-for, not replaced",
    does: "Every artifact carries provenance, uncertainty registers and non-claims so qualified researchers can review the full evidence chain.",
    not: "NEOFORGE does not substitute expert review; review is a human step outside the software." },
  { step: "9 \u00b7 Experimental development & administration", target: "Outside NEOFORGE entirely",
    does: "Nothing \u2014 and the interface says so on every relevant screen.",
    not: "No claim of creating, validating or administering a vaccine for humans or animals." },
];

function viewRosie() {
  app.innerHTML = `
    <h1>Rosie pathway mapping</h1>
    <p class="sub">The reported case of Paul Conyngham and his dog Rosie \u2014 a computationally investigated,
    personalised cancer-vaccine case \u2014 is this project's inspiration. This page maps each reported step of that
    research pathway to what NEOFORGE implements, what it explicitly does not do, and where human expertise and
    the laboratory take over.</p>
    <div class="notice"><span class="ttl">What this page is and is not.</span> The left column is an
    <em>informal narrative mapping</em> based on public reporting of the case, not a verified reproduction of its
    data or methods. NEOFORGE runs its own demonstration cases (a human truth set and an open melanoma case);
    it does not contain Rosie's tumour data, and no outcome \u2014 hers or anyone's \u2014 is claimed, implied or
    predicted by anything in this software.</div>
    <div class="card">
      <table><thead><tr><th>Reported step</th><th>NEOFORGE equivalent</th><th>What NEOFORGE does</th>
        <th>What NEOFORGE does not do</th></tr></thead><tbody>
        ${ROSIE_STEPS.map((s) => `
          <tr><td class="small"><strong>${esc(s.step)}</strong></td>
          <td class="small">${esc(s.target)}</td>
          <td class="small">${esc(s.does)}</td>
          <td class="small dim">${esc(s.not)}</td></tr>`).join("")}
      </tbody></table>
    </div>
    <div class="card">
      <h3 style="margin-top:0">The human-HLA vs canine-DLA gap</h3>
      <p class="small">Dogs have dog leukocyte antigen (DLA) haplotypes, not human HLA. The presentation predictor
      integrated here (MHCflurry) models human class-I alleles only. A Rosie-like canine case therefore cannot be
      run faithfully through this pipeline today: HLA-restricted design on a dog's DLA-typed tumour would be
      unsupported-by-evidence prediction, which the build stage explicitly refuses rather than imputes. Supporting
      canine DLA alleles would require a validated canine predictor \u2014 an honest research gap, recorded here and
      in the construct artifact's uncertainty register rather than papered over.</p>
      <p class="small dim">The reported Rosie pathway also depended on things no software can supply: qualified
      researchers, veterinarians, a research institution's experimental capability, and regulatory oversight.
      NEOFORGE's contribution is confined to the computable prefix of that pathway \u2014 and those boundaries appear
      verbatim in every artifact's non-claims.</p>
    </div>`;
}

/* ---------------- boot ---------------- */

router();
