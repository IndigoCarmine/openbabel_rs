//! Cartesian geometry primitives shared by every force-field energy term.
//!
//! Coordinates are stored flat (`[x0, y0, z0, x1, y1, z1, …]`, one triple per
//! atom in index order) — the same layout OpenBabel uses and that the term
//! export hands us. The `*_at` helpers read one atom's position out of that
//! buffer; the geometric functions return the scalar quantity (bond length,
//! angle, dihedral, out-of-plane angle) an energy term needs.
//!
//! Angles are returned in **radians**. The definitions mirror OpenBabel's own
//! force-field math (`OBForceField::Vector*` in `forcefield.cpp`) so that the
//! energies computed here reproduce OpenBabel's to floating-point tolerance.

/// A 3-vector.
pub(crate) type V3 = [f64; 3];

#[inline]
pub(crate) fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
pub(crate) fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
pub(crate) fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
pub(crate) fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}

/// Read atom `i`'s position from a flat coordinate buffer.
#[inline]
pub(crate) fn pos(coords: &[f64], i: usize) -> V3 {
    [coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]]
}

/// Bond length `|rᵢ − rⱼ|`.
#[inline]
pub(crate) fn distance(coords: &[f64], i: usize, j: usize) -> f64 {
    norm(sub(pos(coords, i), pos(coords, j)))
}

/// Squared bond length (avoids a `sqrt` where only the square is needed).
#[inline]
pub(crate) fn distance_sq(coords: &[f64], i: usize, j: usize) -> f64 {
    let d = sub(pos(coords, i), pos(coords, j));
    dot(d, d)
}

/// Valence angle i–j–k at vertex `j`, in radians (`[0, π]`).
pub(crate) fn angle(coords: &[f64], i: usize, j: usize, k: usize) -> f64 {
    let u = sub(pos(coords, i), pos(coords, j));
    let v = sub(pos(coords, k), pos(coords, j));
    let lu = norm(u);
    let lv = norm(v);
    if lu < 2e-6 || lv < 2e-6 {
        return 0.0;
    }
    let mut cos_t = dot(u, v) / (lu * lv);
    cos_t = cos_t.clamp(-1.0, 1.0);
    cos_t.acos()
}

/// Dihedral angle i–j–k–l in radians (`(-π, π]`).
///
/// Uses the `atan2` formulation (numerically robust near 0/π). Sign follows the
/// right-hand rule about the j–k axis; energy terms depend only on `cos(nφ)`, so
/// the sign convention is immaterial to the energy.
pub(crate) fn dihedral(coords: &[f64], i: usize, j: usize, k: usize, l: usize) -> f64 {
    let b1 = sub(pos(coords, j), pos(coords, i));
    let b2 = sub(pos(coords, k), pos(coords, j));
    let b3 = sub(pos(coords, l), pos(coords, k));
    let n1 = cross(b1, b2);
    let n2 = cross(b2, b3);
    let m = cross(n1, b2);
    let lb2 = norm(b2);
    let y = dot(m, n2) / lb2;
    let x = dot(n1, n2);
    if x.abs() < 1e-12 && y.abs() < 1e-12 {
        return 0.0;
    }
    y.atan2(x)
}

/// Out-of-plane (Wilson) angle of the j→l bond relative to the i–j–k plane,
/// with `j` the central atom, in radians.
///
/// Mirrors the angle OpenBabel's `VectorOOPDerivative` computes: with the three
/// unit bond vectors from the central atom, `sin ψ = (ĵi × ĵk)·ĵl / sin θᵢⱼₖ`.
pub(crate) fn wilson(coords: &[f64], i: usize, j: usize, k: usize, l: usize) -> f64 {
    let c = pos(coords, j);
    let ji = sub(pos(coords, i), c);
    let jk = sub(pos(coords, k), c);
    let jl = sub(pos(coords, l), c);
    let (lji, ljk, ljl) = (norm(ji), norm(jk), norm(jl));
    if lji < 2e-6 || ljk < 2e-6 || ljl < 2e-6 {
        return 0.0;
    }
    let ji = [ji[0] / lji, ji[1] / lji, ji[2] / lji];
    let jk = [jk[0] / ljk, jk[1] / ljk, jk[2] / ljk];
    let jl = [jl[0] / ljl, jl[1] / ljl, jl[2] / ljl];
    let normal = cross(ji, jk);
    let cos_theta = dot(ji, jk).clamp(-1.0, 1.0);
    let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
    if sin_theta < 2e-6 {
        return 0.0;
    }
    (dot(normal, jl) / sin_theta).clamp(-1.0, 1.0).asin()
}

/// Angle between two vectors in **degrees**, matching OpenBabel's `vectorAngle`
/// (including its `±0.9999999` clamp on the cosine).
fn vector_angle_deg(u: V3, v: V3) -> f64 {
    let lu = norm(u);
    let lv = norm(v);
    if lu < 1e-12 || lv < 1e-12 {
        return 0.0;
    }
    let mut dp = dot(u, v) / (lu * lv);
    dp = dp.clamp(-0.999_999_9, 0.999_999_9);
    dp.acos() * (180.0 / std::f64::consts::PI)
}

/// Out-of-plane angle in **degrees**, matching OpenBabel's
/// `Point2PlaneAngle(a, b, c, d)` used by MM2's out-of-plane term. The four atom
/// indices map to OpenBabel's `(a, b, c, d)`: the angle of point `a` above the
/// plane fixed by `b, c, d`, computed as `90° − ∠(normal, a−c)` where
/// `normal = (b−c) × (c−d)`.
pub(crate) fn point2plane_angle(coords: &[f64], a: usize, b: usize, c: usize, d: usize) -> f64 {
    let (pa, pb, pc, pd) = (pos(coords, a), pos(coords, b), pos(coords, c), pos(coords, d));
    let ac = sub(pa, pc);
    let bc = sub(pb, pc);
    let cd = sub(pc, pd);
    let normal = cross(bc, cd);
    90.0 - vector_angle_deg(normal, ac)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn distance_basic() {
        let coords = [0.0, 0.0, 0.0, 3.0, 4.0, 0.0];
        assert!((distance(&coords, 0, 1) - 5.0).abs() < 1e-12);
        assert!((distance_sq(&coords, 0, 1) - 25.0).abs() < 1e-12);
    }

    #[test]
    fn angle_right_and_straight() {
        // 90°: i=(1,0,0) j=(0,0,0) k=(0,1,0)
        let c = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        assert!((angle(&c, 0, 1, 2) - PI / 2.0).abs() < 1e-12);
        // 180°: i=(1,0,0) j=(0,0,0) k=(-1,0,0)
        let c = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0];
        assert!((angle(&c, 0, 1, 2) - PI).abs() < 1e-9);
    }

    #[test]
    fn dihedral_cis_trans() {
        // trans (180°): i and l on opposite sides of the j-k axis
        // j=(0,0,0) k=(1,0,0); i=(0,1,0); l=(1,-1,0)
        let c = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, -1.0, 0.0];
        assert!((dihedral(&c, 0, 1, 2, 3).abs() - PI).abs() < 1e-9);
        // cis (0°): l=(1,1,0)
        let c = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0];
        assert!(dihedral(&c, 0, 1, 2, 3).abs() < 1e-9);
    }

    #[test]
    fn wilson_planar_is_zero() {
        // All four atoms in the z=0 plane → out-of-plane angle 0.
        let c = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, -1.0, 0.0];
        assert!(wilson(&c, 0, 1, 2, 3).abs() < 1e-9);
        // Lift l out of the plane → nonzero.
        let c = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        assert!(wilson(&c, 0, 1, 2, 3).abs() > 0.1);
    }
}
