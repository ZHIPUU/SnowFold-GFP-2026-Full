import React, {useEffect, useMemo, useRef, useState} from 'react'
import {createRoot} from 'react-dom/client'
import Plotly from 'plotly.js-dist-min'
import cytoscape from 'cytoscape'
import fcose from 'cytoscape-fcose'
import {Search, Database, GitBranch, FileText, FlaskConical, Trophy, Network, BarChart3} from 'lucide-react'
import './style.css'

cytoscape.use(fcose)
const API = import.meta.env.VITE_API_URL || 'http://localhost:18000'
const get = (p) => fetch(API + p).then(r => r.json())

function Metric({label, value, tone='gold'}) { return <div className={`metric ${tone}`}><div>{label}</div><b>{value}</b></div> }

function LeftNav({rounds, selectedRound, setSelectedRound, page, setPage}) {
  const items = [
    ['dashboard','总览',Trophy], ['sequences','序列库',Database], ['network','谱系网络',GitBranch], ['plots','图表',BarChart3], ['docs','轮次文档',FileText]
  ]
  return <aside className="left">
    <div className="brand"><FlaskConical size={28}/><div><b>Protein Design Atlas</b><span>GFP evolution cockpit</span></div></div>
    <nav>{items.map(([id,label,Icon]) => <button key={id} onClick={()=>setPage(id)} className={page===id?'active':''}><Icon size={18}/>{label}</button>)}</nav>
    <div className="round-picker"><h3>Rounds</h3><button className={!selectedRound?'selected':''} onClick={()=>setSelectedRound('')}>全部轮次</button>{rounds.map(r=><button key={r.round_key} className={selectedRound===r.round_key?'selected':''} onClick={()=>setSelectedRound(r.round_key)}><span>{r.round_key}</span><em>{r.best_score?Number(r.best_score).toFixed(4):'—'}</em></button>)}</div>
  </aside>
}

function RightPanel({selected, stats}) { return <aside className="right">
  <div className="panel-title">Inspector</div>
  {selected ? <pre>{JSON.stringify(selected,null,2)}</pre> : <div className="empty">点击序列、网络节点或文档条目查看详情。</div>}
  <div className="api-box"><div>API</div><code>GET {API}/api/stats</code><code>GET {API}/api/sequences?min_score=0.94</code><code>GET {API}/api/graph/lineage</code></div>
  {stats?.best && <div className="best-card"><b>Current Best</b><span>{stats.best.best_score?.toFixed?.(4) || stats.best.best_score}</span><small>{stats.best.source_round}</small></div>}
  </aside> }

function Dashboard({stats, trend, scatter}) {
  const trendRef=useRef(); const scatterRef=useRef()
  useEffect(()=>{ if(trend?.length) Plotly.newPlot(trendRef.current,[{x:trend.map(d=>d.round_key), y:trend.map(d=>d.best_score), type:'scatter', mode:'lines+markers', line:{color:'#d4af37', width:4}, marker:{size:10}}], darkLayout('Best score by round')) },[trend])
  useEffect(()=>{ if(scatter?.length) Plotly.newPlot(scatterRef.current,[{x:scatter.map(d=>d.ptm), y:scatter.map(d=>d.chromo), text:scatter.map(d=>`${d.source_round}<br>${d.score?.toFixed?.(4)||d.score}`), mode:'markers', type:'scatter', marker:{size:scatter.map(d=>Math.max(5,(d.score||0.8)*12)), color:scatter.map(d=>d.score), colorscale:'Viridis', showscale:true}}], {...darkLayout('pTM × Chromophore pLDDT'), yaxis:{...darkLayout('').yaxis, range:[0,1], title:'Chromophore pLDDT (0–1)'}, xaxis:{...darkLayout('').xaxis, range:[0,1], title:'pTM'}}) },[scatter])
  return <main className="main"><section className="hero"><span className="eyebrow">GFP Design Lineage</span><h1>从 0.9321 到 0.9449 的实验航迹</h1><p>整合 25 个轮次、二十多万条去重序列、四十万条指标记录，追踪每一次 ProteinMPNN/ESMFold 决策。</p></section><div className="metrics"><Metric label="Sequences" value={stats?.sequences?.toLocaleString?.()}/><Metric label="Metrics" value={stats?.metrics?.toLocaleString?.()} tone="green"/><Metric label="Artifacts" value={stats?.artifacts}/><Metric label="Rounds" value={stats?.rounds} tone="blue"/></div><div className="grid2"><div className="chart" ref={trendRef}/><div className="chart" ref={scatterRef}/></div></main>
}
function darkLayout(title){return {title:{text:title,font:{color:'#f4ecd8'}},paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'rgba(0,0,0,0)',font:{color:'#d9d0bd'},xaxis:{gridcolor:'#2b332f'},yaxis:{gridcolor:'#2b332f'},margin:{l:50,r:20,t:50,b:50}}}

function Sequences({roundKey, onSelect}) { const [rows,setRows]=useState([]); const [q,setQ]=useState('')
  useEffect(()=>{get(`/api/sequences?limit=500${roundKey?'&round_key='+roundKey:''}${q?'&min_score='+q:''}`).then(setRows)},[roundKey,q])
  return <main className="main"><div className="pagehead"><h2>Sequence Vault</h2><div className="search"><Search size={16}/><input placeholder="min score, e.g. 0.94" value={q} onChange={e=>setQ(e.target.value)}/></div></div><table><thead><tr><th>Round</th><th>Score</th><th>pTM</th><th>pLDDT</th><th>Chromo</th><th>Length</th><th>Sequence</th></tr></thead><tbody>{rows.map(r=><tr key={r.id} onClick={()=>onSelect(r)}><td>{r.source_round}</td><td className="score">{r.best_score?.toFixed?.(4)}</td><td>{r.best_ptm?.toFixed?.(4)}</td><td>{r.best_plddt?.toFixed?.(3)}</td><td>{r.best_chromo?.toFixed?.(3)}</td><td>{r.length}</td><td className="seq">{r.sequence}</td></tr>)}</tbody></table></main> }

function NetworkPage({onSelect}) {
  const ref = useRef()
  useEffect(() => {
    let cy
    get('/api/graph/lineage?min_score=0.94&limit=250').then(g => {
      if (!ref.current) return
      cy = cytoscape({
        container: ref.current,
        elements: [...g.nodes, ...g.edges],
        layout: {name:'fcose', animate:false},
        style: [
          {selector:'node', style:{'background-color':'mapData(best_score,0.93,0.946,#315b7c,#e0b84a)','label':'data(label)','color':'#eee','font-size':8,'width':'mapData(best_score,0.93,0.946,18,52)','height':'mapData(best_score,0.93,0.946,18,52)','text-wrap':'wrap'}},
          {selector:'edge', style:{'line-color':'#53645f','target-arrow-color':'#53645f','target-arrow-shape':'triangle','curve-style':'bezier','opacity':0.45}}
        ]
      })
      cy.on('tap','node', e => onSelect(e.target.data()))
    })
    return () => { if (cy) cy.destroy() }
  }, [])
  return <main className="main"><div className="pagehead"><h2>Topology Network</h2><p>父代标签、筛选轮次和高分候选构成的演化网络。</p></div><div className="network" ref={ref}/></main>
}

function DocsPage({roundKey,onSelect}) { const [docs,setDocs]=useState([]); const [body,setBody]=useState(''); useEffect(()=>{get(`/api/documents${roundKey?'?round_key='+roundKey:''}`).then(setDocs)},[roundKey]); return <main className="main docs"><div className="doc-list">{docs.map(d=><button key={d.id} onClick={()=>{onSelect(d); get('/api/documents/'+d.id).then(x=>setBody(x.body))}}>{d.round_key}<b>{d.title}</b><small>{d.path}</small></button>)}</div><article className="markdown"><pre>{body || '选择左侧文档查看完整内容。'}</pre></article></main> }

function PlotsPage(){return <main className="main"><h2>R / Plot Artifacts</h2><p className="muted">R 生成的静态/交互图会写入 artifacts/plots，并由 API 统一索引。当前版本已内置 Plotly 实时图，R worker 可通过 docker compose 生成 SVG/PNG/HTML。</p><div className="plot-note">下一步：运行 <code>docker compose run r-worker Rscript /app/r/scripts/render_all.R</code></div></main>}

function App(){ const [rounds,setRounds]=useState([]), [stats,setStats]=useState(null), [trend,setTrend]=useState([]), [scatter,setScatter]=useState([]), [page,setPage]=useState('dashboard'), [roundKey,setRoundKey]=useState(''), [selected,setSelected]=useState(null)
  useEffect(()=>{get('/api/stats').then(setStats); get('/api/rounds').then(setRounds); get('/api/plots/round-trend').then(setTrend); get('/api/plots/score-scatter?limit=3000').then(setScatter)},[])
  return <div className="app"><LeftNav rounds={rounds} selectedRound={roundKey} setSelectedRound={setRoundKey} page={page} setPage={setPage}/>{page==='dashboard'&&<Dashboard stats={stats} trend={trend} scatter={scatter}/>} {page==='sequences'&&<Sequences roundKey={roundKey} onSelect={setSelected}/>} {page==='network'&&<NetworkPage onSelect={setSelected}/>} {page==='docs'&&<DocsPage roundKey={roundKey} onSelect={setSelected}/>} {page==='plots'&&<PlotsPage/>}<RightPanel selected={selected} stats={stats}/></div> }

createRoot(document.getElementById('root')).render(<App/>);
