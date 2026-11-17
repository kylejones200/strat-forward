//! Athy porosity–depth compaction: φ(z) = φ0 · exp(−z / c).
//!
//! A layer keeps its solid thickness `s`. Buried with its top at depth `z`, its bulk
//! thickness `h` solves  ∫_z^{z+h} (1 − φ(ζ)) dζ = s,  i.e.
//! h − φ0·c·(e^{−z/c} − e^{−(z+h)/c}) = s.
//! Mixed layers use solid-fraction-weighted φ0 and c (an approximation).

/// Porosity at depth `z` below the surface.
pub fn porosity(phi0: f64, c: f64, z: f64) -> f64 {
    phi0 * (-z.max(0.0) / c).exp()
}

/// Bulk thickness of a layer with solid thickness `solid` whose top is at depth `z_top`.
pub fn bulk_thickness(solid: f64, z_top: f64, phi0: f64, c: f64) -> f64 {
    if solid <= 0.0 {
        return 0.0;
    }
    if phi0 <= 0.0 {
        return solid;
    }
    let a = (-z_top.max(0.0) / c).exp();
    // f(h) = h − φ0·c·a·(1 − e^{−h/c}) − s, f'(h) = 1 − φ0·a·e^{−h/c} > 0.
    let mut h = solid / (1.0 - porosity(phi0, c, z_top));
    for _ in 0..50 {
        let e = (-h / c).exp();
        let f = h - phi0 * c * a * (1.0 - e) - solid;
        let df = 1.0 - phi0 * a * e;
        let step = f / df;
        h -= step;
        if step.abs() < 1e-12 * (1.0 + h) {
            break;
        }
    }
    h.max(solid)
}

/// Mean porosity of a layer of bulk thickness `h` and solid thickness `s`.
pub fn mean_porosity(solid: f64, bulk: f64) -> f64 {
    if bulk <= 0.0 {
        0.0
    } else {
        1.0 - solid / bulk
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_thin_surface_layer_has_surface_porosity() {
        let h = bulk_thickness(0.01, 0.0, 0.6, 2000.0);
        assert!((h - 0.01 / 0.4).abs() < 1e-6, "h={h}");
    }

    #[test]
    fn the_solved_thickness_satisfies_the_athy_integral() {
        let (phi0, c) = (0.63, 1960.0);
        for &(s, z) in &[(10.0, 0.0), (50.0, 500.0), (200.0, 2500.0), (1.0, 4000.0)] {
            let h = bulk_thickness(s, z, phi0, c);
            // Midpoint quadrature of ∫(1−φ) over the layer.
            let n = 20_000;
            let dz = h / n as f64;
            let solid: f64 = (0..n)
                .map(|i| (1.0 - porosity(phi0, c, z + (i as f64 + 0.5) * dz)) * dz)
                .sum();
            assert!(
                (solid - s).abs() < 1e-6 * s.max(1.0),
                "s={s} z={z} got {solid}"
            );
        }
    }

    #[test]
    fn burial_only_thins_a_layer() {
        let mut last = f64::INFINITY;
        for z in [0.0, 100.0, 500.0, 1000.0, 3000.0, 6000.0] {
            let h = bulk_thickness(20.0, z, 0.49, 3700.0);
            assert!(h <= last && h >= 20.0);
            last = h;
        }
    }

    #[test]
    fn zero_porosity_is_solid() {
        assert_eq!(bulk_thickness(7.0, 300.0, 0.0, 1000.0), 7.0);
    }
}
