#!/usr/bin/env python3
"""Tiny localhost UI/API bridge for Polyjam. Standard library only."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import secrets
import subprocess
import threading
import time
from http import HTTPStatus
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse, unquote

ROOT = Path(__file__).resolve().parents[1]
WEB = ROOT / "web"
RESULTS = WEB / "results"
BUILTINS = {"tetra", "cube", "octa", "icosa", "dodeca"}
JOBS: dict[str, dict] = {}
JOBS_LOCK = threading.Lock()


def clipped_number(value, lo, hi, kind=float):
    x = kind(value)
    if x < lo or x > hi:
        raise ValueError(f"value {x} outside [{lo}, {hi}]")
    return x


def run_job(job_id: str, config: dict, binary: Path) -> None:
    output_name = f"job-{job_id}.json"
    output_path = RESULTS / output_name
    cmd = [
        str(binary), "search",
        "--piece", config["piece"],
        "--shell", config["shell"],
        "--count", str(config["count"]),
        "--seconds", str(config["seconds"]),
        "--threads", str(config["threads"]),
        "--seed", str(config["seed"]),
        "--tolerance", str(config["tolerance"]),
        "--clearance", str(config["clearance"]),
        "--output", str(output_path),
        "--progress-json",
    ]
    if config.get("startScale", 0) > 0:
        cmd += ["--start-scale", str(config["startScale"])]

    with JOBS_LOCK:
        JOBS[job_id].update(status="running", command=cmd, started=time.time())

    log: list[str] = []
    try:
        proc = subprocess.Popen(
            cmd,
            cwd=ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert proc.stdout is not None
        for raw in proc.stdout:
            line = raw.strip()
            if not line:
                continue
            try:
                event = json.loads(line)
                if isinstance(event, dict) and event.get("type") in {"progress", "done"}:
                    with JOBS_LOCK:
                        JOBS[job_id]["latest"] = event
                    continue
            except json.JSONDecodeError:
                pass
            log.append(line)
            log[:] = log[-20:]
            with JOBS_LOCK:
                JOBS[job_id]["log"] = list(log)
        rc = proc.wait()
        exists = output_path.exists()
        with JOBS_LOCK:
            JOBS[job_id].update(
                status="done" if exists else "failed",
                returnCode=rc,
                result=f"/results/{output_name}" if exists else None,
                finished=time.time(),
                log=list(log),
            )
    except Exception as exc:  # surfaced to the local UI
        with JOBS_LOCK:
            JOBS[job_id].update(status="failed", error=str(exc), finished=time.time(), log=list(log))


class Handler(SimpleHTTPRequestHandler):
    server_version = "Polyjam/1.0"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(WEB), **kwargs)

    def _json(self, payload, status=HTTPStatus.OK):
        body = json.dumps(payload, separators=(",", ":")).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        path = urlparse(self.path).path
        if path == "/api/shapes":
            return self._json({"shapes": sorted(BUILTINS)})
        if path == "/api/results":
            RESULTS.mkdir(parents=True, exist_ok=True)
            items = []
            for p in sorted(RESULTS.glob("*.json"), key=lambda p: p.stat().st_mtime, reverse=True):
                items.append({"name": p.name, "url": f"/results/{p.name}", "mtime": p.stat().st_mtime})
            return self._json({"results": items})
        if path.startswith("/api/jobs/"):
            job_id = unquote(path.removeprefix("/api/jobs/"))
            with JOBS_LOCK:
                job = JOBS.get(job_id)
                if job is None:
                    return self._json({"error": "unknown job"}, HTTPStatus.NOT_FOUND)
                return self._json(job)
        return super().do_GET()

    def do_POST(self):
        path = urlparse(self.path).path
        if path != "/api/search":
            return self._json({"error": "not found"}, HTTPStatus.NOT_FOUND)
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if length <= 0 or length > 65536:
                raise ValueError("invalid request size")
            req = json.loads(self.rfile.read(length))
            piece = str(req.get("piece", "cube"))
            shell = str(req.get("shell", "cube"))
            if piece not in BUILTINS or shell not in BUILTINS:
                raise ValueError("browser search only accepts builtin convex shapes; use the CLI for OBJ files")
            config = {
                "piece": piece,
                "shell": shell,
                "count": clipped_number(req.get("count", 9), 1, 80, int),
                "seconds": clipped_number(req.get("seconds", 10), 0.25, 3600),
                "threads": clipped_number(req.get("threads", max(1, os.cpu_count() or 1)), 1, 128, int),
                "seed": clipped_number(req.get("seed", 1), 0, 2**63 - 1, int),
                "tolerance": clipped_number(req.get("tolerance", 1e-5), 1e-10, 0.1),
                "clearance": clipped_number(req.get("clearance", 0), 0, 1),
                "startScale": clipped_number(req.get("startScale", 0), 0, 1000),
            }
            binary = self.server.binary  # type: ignore[attr-defined]
            if not binary.exists():
                raise ValueError(f"search binary not found at {binary}; build it first")
            job_id = f"{int(time.time())}-{secrets.token_hex(3)}"
            with JOBS_LOCK:
                JOBS[job_id] = {"id": job_id, "status": "queued", "config": config, "latest": None, "log": []}
            threading.Thread(target=run_job, args=(job_id, config, binary), daemon=True).start()
            return self._json({"job": job_id}, HTTPStatus.ACCEPTED)
        except (ValueError, TypeError, json.JSONDecodeError) as exc:
            return self._json({"error": str(exc)}, HTTPStatus.BAD_REQUEST)

    def log_message(self, fmt, *args):
        print(f"[{self.log_date_time_string()}] {fmt % args}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8080)
    ap.add_argument("--binary", type=Path, default=ROOT / "build" / "polyjam")
    args = ap.parse_args()
    RESULTS.mkdir(parents=True, exist_ok=True)
    server = ThreadingHTTPServer((args.host, args.port), Handler)
    server.binary = args.binary.resolve()  # type: ignore[attr-defined]
    print(f"Polyjam UI: http://{args.host}:{args.port}/")
    print(f"Search binary: {server.binary}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
