"""Deterministic rigid start placement, before unconstrained surface relaxation.

This changes neither molecular topology nor the race force model. Ambiguous
principal axes do not define a vehicle front. The small Jacobi solve is also
implemented in the browser C++ path, with identical thresholds and tie rules.
"""
from __future__ import annotations
import math
import numpy as np

POSE_VERSION = "rigid-start-pose-v1"
DEFAULT_START_POSE = {"mode": "auto", "yaw_deg": 0., "pitch_deg": 0., "roll_deg": 0., "axis_atoms": []}
PRESERVE_START_POSE = {**DEFAULT_START_POSE, "mode": "preserve"}
HEADING_RATIO = 1.3
PLANE_RATIO = .25


def validate_start_pose(value=None, *, atom_count=None):
    if value is None:
        value = DEFAULT_START_POSE
    if not isinstance(value, dict) or set(value) - set(DEFAULT_START_POSE):
        raise ValueError("起跑姿态参数格式不正确。")
    result = {**DEFAULT_START_POSE, **value}
    if result["mode"] not in ("auto", "manual", "preserve"):
        raise ValueError("起跑摆放方式必须是 auto、manual 或 preserve。")
    for key in ("yaw_deg", "pitch_deg", "roll_deg"):
        v = result[key]
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or not -180 <= v <= 180:
            raise ValueError("起跑旋转角度必须是 −180° 到 180° 的有限数值。")
        result[key] = float(v)
    ids = result["axis_atoms"]
    if (not isinstance(ids, list) or len(ids) not in (0, 2)
            or any(isinstance(v, bool) or not isinstance(v, int) or not 1 <= v <= 256 for v in ids)
            or (len(ids) == 2 and ids[0] == ids[1])):
        raise ValueError("车头参考需要两个不同的有效原子编号，按车尾、车头顺序填写。")
    if atom_count is not None and any(v > atom_count for v in ids):
        raise ValueError("车头参考原子不属于当前分子，请重新选择。")
    result["axis_atoms"] = list(ids)
    # Normalize unused controls so the requested pose has unambiguous semantics.
    if result["mode"] != "manual":
        for key in ("yaw_deg", "pitch_deg", "roll_deg"):
            result[key] = 0.
    if result["mode"] != "auto":
        result["axis_atoms"] = []
    return result


def _principal_axes(x):
    center = [sum(float(row[k]) for row in x) / len(x) for k in range(3)]
    covariance = [[sum((float(row[i]) - center[i]) * (float(row[j]) - center[j]) for row in x) / len(x)
                   for j in range(3)] for i in range(3)]
    vectors = [[float(i == j) for j in range(3)] for i in range(3)]
    for _ in range(32):
        p, q = max(((0, 1), (0, 2), (1, 2)), key=lambda pair: abs(covariance[pair[0]][pair[1]]))
        scale = max(1., sum(abs(covariance[i][i]) for i in range(3)))
        if abs(covariance[p][q]) <= 1e-14 * scale:
            break
        angle = .5 * math.atan2(2 * covariance[p][q], covariance[q][q] - covariance[p][p])
        c, s = math.cos(angle), math.sin(angle)
        pp, qq, pq = covariance[p][p], covariance[q][q], covariance[p][q]
        for k in range(3):
            if k not in (p, q):
                kp, kq = covariance[k][p], covariance[k][q]
                covariance[k][p] = covariance[p][k] = c * kp - s * kq
                covariance[k][q] = covariance[q][k] = s * kp + c * kq
            vp, vq = vectors[k][p], vectors[k][q]
            vectors[k][p], vectors[k][q] = c * vp - s * vq, s * vp + c * vq
        covariance[p][p] = c*c*pp - 2*s*c*pq + s*s*qq
        covariance[q][q] = s*s*pp + 2*s*c*pq + c*c*qq
        covariance[p][q] = covariance[q][p] = 0.
    order = sorted(range(3), key=lambda i: -covariance[i][i])
    return ([max(0., covariance[i][i]) for i in order],
            [np.array([vectors[k][i] for k in range(3)]) for i in order])


def _minimum_rotation(a, b):
    a, b = a / np.linalg.norm(a), b / np.linalg.norm(b)
    cosine = float(np.dot(a, b))
    if cosine > 1 - 1e-12:
        return np.eye(3)
    if cosine < -1 + 1e-12:
        index = int(np.argmin(np.abs(a)))
        seed = np.eye(3)[index]
        axis = np.cross(a, seed); axis /= np.linalg.norm(axis)
        return 2 * np.outer(axis, axis) - np.eye(3)
    axis = np.cross(a, b)
    matrix = np.array([[0., -axis[2], axis[1]], [axis[2], 0., -axis[0]], [-axis[1], axis[0], 0.]])
    return np.eye(3) + matrix + matrix @ matrix / (1 + cosine)


def place_start_pose(positions, masses, elements, value=None, *, minimum_z=2.9):
    x = np.asarray(positions, dtype=np.float64)
    masses = np.asarray(masses, dtype=np.float64)
    config = validate_start_pose(value, atom_count=len(x))
    center = np.average(x, axis=0, weights=masses)
    rotation = np.eye(3)
    reference = [i for i, element in enumerate(elements) if element != "H"] or list(range(len(x)))
    points = x[reference]
    variances, axes = _principal_axes(points)
    heading = variances[0] > 1e-12 and variances[0] >= HEADING_RATIO * variances[1]
    plane = variances[1] > 1e-12 and variances[2] <= PLANE_RATIO * variances[1]
    strategy = "input-preserved"
    heading_ids = []
    if config["mode"] == "manual":
        yaw, pitch, roll = [math.radians(config[k]) for k in ("yaw_deg", "pitch_deg", "roll_deg")]
        cy, sy, cp, sp, cr, sr = math.cos(yaw), math.sin(yaw), math.cos(pitch), math.sin(pitch), math.cos(roll), math.sin(roll)
        rotation = np.array([[cy*cp, cy*sp*sr-sy*cr, cy*sp*cr+sy*sr],
                             [sy*cp, sy*sp*sr+cy*cr, sy*sp*cr-cy*sr], [-sp, cp*sr, cp*cr]])
        strategy = "manual-file-axes-roll-pitch-yaw"
    elif config["mode"] == "auto":
        forward, normal = axes[0].copy(), axes[2].copy()
        if normal[2] < -1e-12 or (abs(normal[2]) <= 1e-12 and next((v for v in normal if abs(v) > 1e-12), 1.) < 0):
            normal *= -1
        ids = config["axis_atoms"]
        if ids:
            forward = x[ids[1] - 1] - x[ids[0] - 1]
            if np.linalg.norm(forward) <= 1e-8:
                raise ValueError("车头参考原子距离过近，无法定义方向。")
            forward /= np.linalg.norm(forward)
            heading = True; heading_ids = ids
        elif heading:
            projections = (points - points.mean(axis=0)) @ forward
            # A principal axis has no physical front/back polarity. Pick its
            # closest +x direction, rather than turn a nearly aligned car 180°
            # because of arbitrary eigensolver or atom-order signs.
            sign = next((v for v in forward if abs(v) > 1e-12), 1.)
            if sign < 0:
                forward *= -1; projections *= -1
            low, high = float(projections.min()), float(projections.max())
            tolerance = 1e-8 * max(1., high - low)
            rear = next(i for i, p in enumerate(projections) if p <= low + tolerance)
            front = next(i for i, p in enumerate(projections) if p >= high - tolerance)
            heading_ids = [reference[rear] + 1, reference[front] + 1]
        if heading and plane and np.linalg.norm(np.cross(normal, forward)) > 1e-8:
            right = np.cross(normal, forward); right /= np.linalg.norm(right)
            up = np.cross(forward, right)
            rotation = np.array([forward, right, up])
            strategy = "automatic-plane-and-heading"
        elif heading:
            rotation = _minimum_rotation(forward, np.array([1., 0., 0.]))
            strategy = "automatic-heading-preserve-roll"
        elif plane:
            rotation = _minimum_rotation(normal, np.array([0., 0., 1.]))
            strategy = "automatic-plane-preserve-yaw"
        else:
            strategy = "ambiguous-geometry-input-preserved"
    # Avoid altering the legacy preserve path by subtracting/readding COM.
    result = (np.array(x, copy=True) if np.array_equal(rotation, np.eye(3))
              else (x - center) @ rotation.T + center)
    result[:, 0] -= center[0]; result[:, 1] -= center[1]
    result[:, 2] += minimum_z - result[:, 2].min()
    report = {"protocol_version": POSE_VERSION, "requested": config, "strategy": strategy,
              "rotation_matrix": rotation.tolist(), "principal_variances_A2": variances,
              "heading_confident": bool(heading), "plane_confident": bool(plane),
              "heading_atom_ids": heading_ids,
              "qualification": "Rigid pre-relaxation placement only; unconstrained surface relaxation and subsequent dynamics may change orientation."}
    return result, report


def relaxation_heading_report(positions, report):
    ids = report.get("heading_atom_ids", [])
    if not ids:
        return {"heading_defined": False}
    direction = positions[ids[1] - 1] - positions[ids[0] - 1]
    return {"heading_defined": True, "heading_atom_ids": ids,
            "yaw_deg": math.degrees(math.atan2(direction[1], direction[0])),
            "tilt_deg": math.degrees(math.atan2(direction[2], math.hypot(direction[0], direction[1])))}
