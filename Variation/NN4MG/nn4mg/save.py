from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Any, Mapping

import torch


def _jsonify(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(k): _jsonify(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonify(v) for v in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if isinstance(value, Path):
        return str(value)
    return str(value)


def prepare_run_dir(base_dir: str | Path, run_name: str | None = None) -> Path:
    base_path = Path(base_dir).expanduser()
    base_path.mkdir(parents=True, exist_ok=True)
    if run_name is None:
        run_name = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = base_path / run_name
    if run_dir.exists():
        suffix = 1
        while True:
            candidate = base_path / f"{run_name}_{suffix:02d}"
            if not candidate.exists():
                run_dir = candidate
                break
            suffix += 1
    run_dir.mkdir(parents=True, exist_ok=False)
    return run_dir


def save_run(
    run_dir: str | Path,
    net: torch.nn.Module,
    config: Mapping[str, Any],
    model_filename: str = "model.pt",
    config_filename: str = "run_config.json",
) -> Path:
    run_path = Path(run_dir).expanduser()
    run_path.mkdir(parents=True, exist_ok=True)
    torch.save(net.state_dict(), run_path / model_filename)
    config_out = _jsonify(config)
    with (run_path / config_filename).open("w", encoding="utf-8") as f:
        json.dump(config_out, f, indent=2, sort_keys=True, ensure_ascii=True)
    return run_path
