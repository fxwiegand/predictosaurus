use colored::Colorize;
use std::fmt;

#[derive(Debug, PartialEq, Eq, Ord, PartialOrd, Clone, Copy)]
pub(crate) enum Impact {
    None,
    Modifier,
    Low,
    Moderate,
    High,
}

impl fmt::Display for Impact {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let colored_str = match *self {
            Impact::None => "None".green(),
            Impact::Modifier => "Modifier".green(),
            Impact::Low => "Low".yellow(),
            Impact::Moderate => "Moderate".red(),
            Impact::High => "High".purple(),
        };
        write!(f, "{}", colored_str)
    }
}

#[cfg(test)]
mod impact_enum_tests {
    use super::Impact;
    use colored::Colorize;

    #[test]
    fn test_impact_order() {
        assert!(Impact::High > Impact::Moderate);
        assert!(Impact::High > Impact::Low);
        assert!(Impact::High > Impact::Modifier);
        assert!(Impact::High > Impact::None);
        assert!(Impact::Moderate > Impact::Low);
        assert!(Impact::Moderate > Impact::Modifier);
        assert!(Impact::Moderate > Impact::None);
        assert!(Impact::Low > Impact::Modifier);
        assert!(Impact::Low > Impact::None);
        assert!(Impact::Modifier > Impact::None);
    }

    #[test]
    fn display_none_impact() {
        let impact = Impact::None;
        assert_eq!(format!("{}", impact), "None".green().to_string());
    }

    #[test]
    fn display_modifier_impact() {
        let impact = Impact::Modifier;
        assert_eq!(format!("{}", impact), "Modifier".green().to_string());
    }

    #[test]
    fn display_low_impact() {
        let impact = Impact::Low;
        assert_eq!(format!("{}", impact), "Low".yellow().to_string());
    }

    #[test]
    fn display_moderate_impact() {
        let impact = Impact::Moderate;
        assert_eq!(format!("{}", impact), "Moderate".red().to_string());
    }

    #[test]
    fn display_high_impact() {
        let impact = Impact::High;
        assert_eq!(format!("{}", impact), "High".purple().to_string());
    }
}
