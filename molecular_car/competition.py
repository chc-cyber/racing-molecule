"""Versioned race rules, fixed before molecular/surface preparation."""
import copy
import math
from .start_pose import validate_start_pose

DEFAULT_COMPETITION = dict(mode="standard", version="race-v1", length_nm=20,
    temperature=150, strength=.6, surface="fourier-au", lattice_scale=1.5,
    corrugation=.08, growth=1, timestep_fs=1, damping_fs=350, random_seed=98245165,
    start_pose=validate_start_pose())


def validate_competition(data=None):
    if data is None:
        data = {}
    if not isinstance(data, dict) or set(data) - set(DEFAULT_COMPETITION):
        raise ValueError("赛事配置包含未知参数。")
    result = copy.deepcopy(DEFAULT_COMPETITION)
    result.update(copy.deepcopy(data))
    if result["mode"] not in ("standard", "custom"):
        raise ValueError("请选择标准赛事或自定义赛事。")
    for key in ("version", "length_nm", "timestep_fs", "damping_fs", "random_seed"):
        if type(result[key]) is not type(DEFAULT_COMPETITION[key]) or result[key] != DEFAULT_COMPETITION[key]:
            raise ValueError("当前赛事版本固定为 20 nm、1 fs 步长与 350 fs 阻尼。")
    if result["surface"] not in ("fourier-au", "periodic-au", "continuum"):
        raise ValueError("未知表面模型。")
    for key, low, high in (("temperature", 0, 1000), ("strength", 0, 1),
                           ("lattice_scale", 0, 3), ("corrugation", 0, .6), ("growth", .1, 2)):
        value = result[key]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not low <= value <= high:
            raise ValueError(f"{key} 必须在 {low}–{high} 范围内。")
    result["start_pose"] = validate_start_pose(result["start_pose"])
    if result["mode"] == "standard" and result != DEFAULT_COMPETITION:
        raise ValueError("标准赛事参数固定；请切换到自定义赛事。")
    return result


def engine_parameters(competition):
    return ({key: competition[key] for key in
             ("temperature", "damping_fs", "corrugation", "lattice_scale", "surface")},
            dict(mode="gradient", strength=competition["strength"], growth=competition["growth"]))
