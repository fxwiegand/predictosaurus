use anyhow::{anyhow, Result};
use colored::Colorize;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::str::FromStr;

#[derive(Debug, PartialEq, Eq, Ord, PartialOrd, Clone, Copy, Serialize, Deserialize)]
pub(crate) enum Impact {
    None,
    Modifier,
    Low,
    Moderate,
    High,
}

impl FromStr for Impact {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        match s.to_lowercase().as_str() {
            "none" => Ok(Impact::None),
            "modifier" => Ok(Impact::Modifier),
            "low" => Ok(Impact::Low),
            "moderate" => Ok(Impact::Moderate),
            "high" => Ok(Impact::High),
            _ => Err(anyhow!("'{}' is not a valid Impact variant", s)),
        }
    }
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
        write!(f, "{colored_str}")
    }
}

impl Impact {
    /// Serializes the Impact variant to a string without color
    pub(crate) fn to_raw_string(self) -> &'static str {
        match self {
            Impact::None => "None",
            Impact::Modifier => "Modifier",
            Impact::Low => "Low",
            Impact::Moderate => "Moderate",
            Impact::High => "High",
        }
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

    #[test]
    fn from_str_none() {
        let impact = "none".parse::<Impact>().unwrap();
        assert_eq!(impact, Impact::None);
    }

    #[test]
    fn from_str_modifier() {
        let impact = "modifier".parse::<Impact>().unwrap();
        assert_eq!(impact, Impact::Modifier);
    }

    #[test]
    fn from_str_low() {
        let impact = "low".parse::<Impact>().unwrap();
        assert_eq!(impact, Impact::Low);
    }

    #[test]
    fn from_str_moderate() {
        let impact = "moderate".parse::<Impact>().unwrap();
        assert_eq!(impact, Impact::Moderate);
    }

    #[test]
    fn from_str_high() {
        let impact = "high".parse::<Impact>().unwrap();
        assert_eq!(impact, Impact::High);
    }

    #[test]
    fn from_str_invalid() {
        let impact = "invalid".parse::<Impact>();
        assert!(impact.is_err());
    }
}
