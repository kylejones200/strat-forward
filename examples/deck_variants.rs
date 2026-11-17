// Reservoir decks from one run of the reference configuration with different export options, to show which
// export choices change the grid. Writes <dir>/<name>.DATA (Eclipse format) and prints each GridSummary.
//
//   cargo run --release --example deck_variants -- out/deck_variants
use strat_forward::reservoir::{self, ReservoirOptions};
use strat_forward::{Config, Model};

fn main() {
    let dir = std::path::PathBuf::from(
        std::env::args()
            .nth(1)
            .unwrap_or("out/deck_variants".into()),
    );
    std::fs::create_dir_all(&dir).unwrap();
    let mut m = Model::new(Config::example()).unwrap();
    m.run().unwrap();
    let base = ReservoirOptions::default();
    let variants: Vec<(&str, ReservoirOptions)> = vec![
        ("default", base.clone()),
        (
            "no_reburial",
            ReservoirOptions {
                overburden_m: Some(0.0),
                ..base.clone()
            },
        ),
        (
            // Same window as the default deck, so only burial differs.
            "no_reburial_window_55km",
            ReservoirOptions {
                overburden_m: Some(0.0),
                x_start_m: Some(55_000.0),
                ..base.clone()
            },
        ),
        (
            "min_thickness_0p5",
            ReservoirOptions {
                min_thickness_m: 0.5,
                ..base.clone()
            },
        ),
        (
            "kv_kh_1",
            ReservoirOptions {
                kv_kh: 1.0,
                ..base.clone()
            },
        ),
        (
            "window_40km",
            ReservoirOptions {
                x_start_m: Some(40_000.0),
                ..base.clone()
            },
        ),
    ];
    for (name, opt) in variants {
        match reservoir::deck(&m, &opt) {
            Ok(d) => {
                std::fs::write(dir.join(format!("{name}.DATA")), &d.text).unwrap();
                println!("{name} {}", serde_json::to_string(&d.summary).unwrap());
            }
            Err(e) => println!("{name} error: {e}"),
        }
    }
}
