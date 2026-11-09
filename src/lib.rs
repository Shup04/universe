//! Experimental wave-field simulation, independent of a window or random heating.
mod simulation;
pub use simulation::{Universe, View};

#[cfg(feature = "gui")]
pub mod gui;
