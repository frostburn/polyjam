const $ = id => document.getElementById(id)
const canvas = $('view')
const ctx = canvas.getContext('2d')

const state = {
  result: null,
  yaw: 0.72,
  pitch: -0.48,
  distance: 7,
  selected: -1,
  drag: null,
  projectedCenters: [],
  currentJob: null,
}

function qrotate(q, v) {
  const [w,x,y,z] = q
  const [vx,vy,vz] = v
  const tx = 2 * (y*vz - z*vy)
  const ty = 2 * (z*vx - x*vz)
  const tz = 2 * (x*vy - y*vx)
  return [
    vx + w*tx + (y*tz-z*ty),
    vy + w*ty + (z*tx-x*tz),
    vz + w*tz + (x*ty-y*tx),
  ]
}
function add(a,b) { return [a[0]+b[0],a[1]+b[1],a[2]+b[2]] }
function scale(v,s) { return [v[0]*s,v[1]*s,v[2]*s] }
function len(v) { return Math.hypot(v[0],v[1],v[2]) }
function cameraPoint(v) {
  const cy=Math.cos(state.yaw), sy=Math.sin(state.yaw)
  const cp=Math.cos(state.pitch), sp=Math.sin(state.pitch)
  const x1=cy*v[0]-sy*v[1]
  const y1=sy*v[0]+cy*v[1]
  const z1=v[2]
  return [x1, cp*y1-sp*z1, sp*y1+cp*z1]
}
function project(v, w, h) {
  const c=cameraPoint(v)
  const z=Math.max(.08, state.distance-c[2])
  const f=Math.min(w,h)*0.72
  return [w/2+c[0]*f/z, h/2-c[1]*f/z, c[2], z]
}

function colorFor(i, alpha=1) {
  const hue=(i*137.508+205)%360
  return `hsla(${hue},68%,62%,${alpha})`
}

function geometryObjects() {
  const r=state.result
  if (!r) return []
  const explode=Number($('explode').value)
  const slice=Number($('slice').value)
  const out=[]
  const shellVerts=r.shell.vertices.map(v => scale(v,r.scale))
  out.push({kind:'shell', index:-1, verts:shellVerts, faces:r.shell.faces})
  state.projectedCenters=[]
  r.poses.forEach((pose,i) => {
    const limit = slice * r.scale
    if (pose.p[2] > limit) return
    const radial = len(pose.p) > 1e-12 ? scale(pose.p, explode*.35) : [0,0,0]
    const center=add(pose.p,radial)
    const verts=r.piece.vertices.map(v => add(qrotate(pose.q,v),center))
    out.push({kind:'piece', index:i, verts, faces:r.piece.faces, center})
  })
  return out
}

function draw() {
  // Keep the render loop alive while result loading is asynchronous. The first
  // frame usually runs before a bundled, server, or local result has loaded.
  requestAnimationFrame(draw)
  const rect=canvas.getBoundingClientRect()
  const dpr=Math.max(1,window.devicePixelRatio||1)
  const W=Math.round(rect.width*dpr), H=Math.round(rect.height*dpr)
  if (canvas.width!==W || canvas.height!==H) { canvas.width=W; canvas.height=H }
  ctx.setTransform(dpr,0,0,dpr,0,0)
  const w=rect.width,h=rect.height
  ctx.clearRect(0,0,w,h)
  if (!state.result) {
    ctx.fillStyle='#708096'; ctx.font='14px system-ui'; ctx.textAlign='center'
    ctx.fillText('Run a search or open a result JSON',w/2,h/2)
    return
  }

  const faces=[]
  const edgeSegments=[]
  const shellAlpha=Number($('shellOpacity').value)
  const showFaces=$('faces').checked, showEdges=$('edges').checked
  const objs=geometryObjects()

  for (const obj of objs) {
    const projected=obj.verts.map(v=>project(v,w,h))
    if (obj.kind==='piece') {
      const p=project(obj.center,w,h)
      state.projectedCenters[obj.index]=p
    }
    for (const face of obj.faces) {
      const pts=face.map(i=>projected[i])
      const depth=pts.reduce((s,p)=>s+p[2],0)/pts.length
      faces.push({obj,pts,depth})
      for (let j=0;j<pts.length;j++) edgeSegments.push({obj,a:pts[j],b:pts[(j+1)%pts.length],depth})
    }
  }
  // Painter's algorithm: far faces first in camera's rotated z coordinate.
  faces.sort((a,b)=>a.depth-b.depth)
  if (showFaces) for (const f of faces) {
    ctx.beginPath(); ctx.moveTo(f.pts[0][0],f.pts[0][1])
    for (let i=1;i<f.pts.length;i++) ctx.lineTo(f.pts[i][0],f.pts[i][1])
    ctx.closePath()
    if (f.obj.kind==='shell') {
      ctx.fillStyle=`rgba(174,194,218,${shellAlpha})`
    } else {
      const selected=f.obj.index===state.selected
      ctx.fillStyle=colorFor(f.obj.index,selected?.70:.38)
    }
    ctx.fill()
  }
  if (showEdges) {
    // De-duplicated-looking edges are not necessary; repeated strokes make contacts pleasantly crisp.
    for (const e of edgeSegments) {
      ctx.beginPath(); ctx.moveTo(e.a[0],e.a[1]); ctx.lineTo(e.b[0],e.b[1])
      if (e.obj.kind==='shell') { ctx.strokeStyle='rgba(205,220,238,.52)'; ctx.lineWidth=1.2 }
      else { ctx.strokeStyle=e.obj.index===state.selected?'rgba(255,255,255,.96)':'rgba(230,239,249,.58)'; ctx.lineWidth=e.obj.index===state.selected?2.2:1.0 }
      ctx.stroke()
    }
  }
}
requestAnimationFrame(draw)

function setResult(r, label='result') {
  if (!r || r.format!=='polyjam-result-v1') throw new Error('Not a polyjam-result-v1 file')
  state.result=r; state.selected=-1
  state.distance=Math.max(3.4, r.scale*2.35)
  $('statScale').textContent=Number(r.scale).toFixed(8)
  $('statCount').textContent=r.count
  $('statFeasible').textContent=r.feasible?'yes':'NO'
  $('statFeasible').className=r.feasible?'good':'bad'
  $('statViolation').textContent=Number(r.metrics?.maxViolation ?? NaN).toExponential(3)
  $('note').textContent=r.note||''
  $('badge').textContent=`${label} · ${r.count} × ${r.piece.name} in ${r.shell.name}`
  $('selected').textContent='Drag to orbit · wheel to zoom · click a piece'
}

async function loadResult(url,label) {
  if (url.startsWith('sample:')) {
    const i=Number(url.slice(7)); const sample=(window.POLYJAM_SAMPLES||[])[i]
    if (!sample) throw new Error('Unknown bundled sample')
    setResult(sample.data,sample.name); return
  }
  const res=await fetch(url,{cache:'no-store'})
  if (!res.ok) throw new Error(`HTTP ${res.status}`)
  setResult(await res.json(),label)
}

async function refreshResults(selectNewest=false) {
  const sel=$('results')
  try {
    const data=await (await fetch('/api/results',{cache:'no-store'})).json()
    const old=sel.value
    sel.innerHTML=''
    for (const item of data.results) {
      const o=document.createElement('option'); o.value=item.url; o.textContent=item.name; sel.appendChild(o)
    }
    if (!data.results.length) { const o=document.createElement('option'); o.textContent='No saved results'; o.value=''; sel.appendChild(o); return }
    if (!selectNewest && [...sel.options].some(o=>o.value===old)) sel.value=old
    else sel.selectedIndex=0
    if (selectNewest || !state.result) await loadResult(sel.value,sel.options[sel.selectedIndex].textContent)
  } catch (e) {
    const samples=window.POLYJAM_SAMPLES||[]
    sel.innerHTML=''
    samples.forEach((sample,i)=>{ const o=document.createElement('option');o.value=`sample:${i}`;o.textContent=`${sample.name} (bundled)`;sel.appendChild(o) })
    if (samples.length && !state.result) await loadResult('sample:0',samples[0].name)
    if (!samples.length) sel.innerHTML='<option>Server unavailable — use Open JSON</option>'
  }
}

async function initShapes() {
  let shapes=['cube','tetra','octa','icosa','dodeca']
  try { shapes=(await (await fetch('/api/shapes')).json()).shapes } catch {}
  for (const id of ['piece','shell']) {
    const s=$(id); s.innerHTML=''
    for (const name of shapes) { const o=document.createElement('option');o.value=o.textContent=name;s.appendChild(o) }
    s.value='cube'
  }
  $('threads').value=Math.min(8,navigator.hardwareConcurrency||4)
}

async function startSearch() {
  const button=$('search'), status=$('job')
  button.disabled=true
  const body={
    piece:$('piece').value,shell:$('shell').value,count:Number($('count').value),seconds:Number($('seconds').value),
    threads:Number($('threads').value),seed:Number($('seed').value),tolerance:Number($('tolerance').value),
    clearance:Number($('clearance').value),startScale:Number($('startScale').value)
  }
  try {
    const res=await fetch('/api/search',{method:'POST',headers:{'content-type':'application/json'},body:JSON.stringify(body)})
    const data=await res.json(); if (!res.ok) throw new Error(data.error||`HTTP ${res.status}`)
    state.currentJob=data.job; status.textContent=`job ${data.job} · queued`
    await pollJob(data.job)
  } catch(e) {
    status.textContent=`error: ${e.message}`; button.disabled=false
  }
}

async function pollJob(id) {
  const button=$('search'), status=$('job')
  while (state.currentJob===id) {
    const res=await fetch(`/api/jobs/${encodeURIComponent(id)}`,{cache:'no-store'})
    const job=await res.json()
    const p=job.latest
    if (p) status.textContent=`${job.status} · scale ${Number(p.scale).toFixed(7)} · maxV ${Number(p.maxViolation).toExponential(2)}`
    else status.textContent=`job ${id} · ${job.status}`
    if (job.status==='done') {
      state.currentJob=null; button.disabled=false
      await refreshResults(true)
      return
    }
    if (job.status==='failed') {
      state.currentJob=null; button.disabled=false
      status.textContent=`failed: ${job.error||job.log?.slice(-1)[0]||'search exited without a result'}`
      return
    }
    await new Promise(r=>setTimeout(r,450))
  }
}

$('search').addEventListener('click',startSearch)
$('refresh').addEventListener('click',()=>refreshResults(false))
$('results').addEventListener('change',e=>e.target.value&&loadResult(e.target.value,e.target.options[e.target.selectedIndex].textContent).catch(err=>alert(err.message)))
$('file').addEventListener('change',async e=>{
  const f=e.target.files[0]; if (!f) return
  try { setResult(JSON.parse(await f.text()),f.name) } catch(err) { alert(err.message) }
})

for (const id of ['shellOpacity','explode','slice']) $(id).addEventListener('input',()=>{
  $('shellOpacityOut').textContent=`${Math.round(Number($('shellOpacity').value)*100)}%`
  $('explodeOut').textContent=`${Math.round(Number($('explode').value)*100)}%`
  $('sliceOut').textContent=`${Math.round((Number($('slice').value)+1)*50)}%`
})

canvas.addEventListener('pointerdown',e=>{
  state.drag={x:e.clientX,y:e.clientY,yaw:state.yaw,pitch:state.pitch,moved:false}; canvas.setPointerCapture(e.pointerId); canvas.classList.add('dragging')
})
canvas.addEventListener('pointermove',e=>{
  if (!state.drag) return
  const dx=e.clientX-state.drag.x,dy=e.clientY-state.drag.y
  if (Math.hypot(dx,dy)>3) state.drag.moved=true
  state.yaw=state.drag.yaw+dx*.008
  state.pitch=Math.max(-1.48,Math.min(1.48,state.drag.pitch+dy*.008))
})
canvas.addEventListener('pointerup',e=>{
  if (state.drag && !state.drag.moved && state.result) {
    const rect=canvas.getBoundingClientRect(); const x=e.clientX-rect.left,y=e.clientY-rect.top
    let best=-1,bestD=28
    state.projectedCenters.forEach((p,i)=>{ if(!p)return; const d=Math.hypot(p[0]-x,p[1]-y);if(d<bestD){bestD=d;best=i} })
    state.selected=best
    if (best>=0) {
      const p=state.result.poses[best]
      $('selected').textContent=`piece ${best+1} · p=[${p.p.map(x=>Number(x).toFixed(4)).join(', ')}]`
    } else $('selected').textContent='No piece selected'
  }
  state.drag=null; canvas.classList.remove('dragging')
})
canvas.addEventListener('wheel',e=>{ e.preventDefault(); state.distance=Math.max(.8,state.distance*Math.exp(e.deltaY*.001)); },{passive:false})
canvas.addEventListener('dblclick',()=>{ state.yaw=.72;state.pitch=-.48;state.distance=Math.max(3.4,(state.result?.scale||3)*2.35) })

Promise.all([initShapes(),refreshResults(true)])
