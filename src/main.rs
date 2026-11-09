fn main() {
    if let Err(error) = universe::gui::run() {
        eprintln!("Unable to run graphical demo: {error:?}");
        std::process::exit(1);
    }
}
