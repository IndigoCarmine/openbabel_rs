//! Chemical reactions: reactant / product / agent molecule lists.

use cxx::UniquePtr;
use openbabel_sys::ffi;

use crate::error::Error;
use crate::mol::Molecule;

/// A chemical reaction — ordered lists of reactant, product, and agent
/// molecules, plus a title, comment, and reversibility flag.
///
/// Read and written through OpenBabel's reaction formats, chiefly `"rsmi"`
/// (reaction SMILES like `"C=C.O>>CCO"`) and `"rxn"` (MDL RXN). For *applying* a
/// reaction as a graph edit, see [`Transform`](crate::Transform) instead; this
/// type models the reaction as data.
///
/// ```no_run
/// use openbabel::Reaction;
/// let rxn = Reaction::parse("C=C.O>>CCO", "rsmi").unwrap();
/// assert_eq!(rxn.num_reactants(), 2);
/// assert_eq!(rxn.num_products(), 1);
/// ```
pub struct Reaction {
    inner: UniquePtr<ffi::Reaction>,
}

impl Reaction {
    /// An empty reaction (no reactants, products, or agents).
    pub fn new() -> Self {
        Reaction {
            inner: crate::with_ob(ffi::reaction_new),
        }
    }

    /// Parse a reaction document in `format` (e.g. `"rsmi"`, `"rxn"`).
    ///
    /// Returns [`Error::Parse`] if `format` is unknown or the text cannot be
    /// parsed as a reaction.
    pub fn parse(text: &str, format: &str) -> Result<Self, Error> {
        let inner = crate::with_ob(|| ffi::reaction_read(format, text));
        if inner.is_null() {
            Err(Error::Parse {
                format: format.to_string(),
            })
        } else {
            Ok(Reaction { inner })
        }
    }

    /// Serialize this reaction in `format` (e.g. `"rsmi"`, `"rxn"`).
    ///
    /// Returns [`Error::UnknownFormat`] if `format` is not a writable format.
    pub fn write(&self, format: &str) -> Result<String, Error> {
        let mut ok = false;
        let out = crate::with_ob(|| ffi::reaction_write(self.as_inner(), format, &mut ok));
        if ok {
            Ok(out)
        } else {
            Err(Error::UnknownFormat {
                format: format.to_string(),
            })
        }
    }

    /// Number of reactant molecules.
    pub fn num_reactants(&self) -> u32 {
        crate::with_ob(|| ffi::reaction_num_reactants(self.as_inner()))
    }

    /// Number of product molecules.
    pub fn num_products(&self) -> u32 {
        crate::with_ob(|| ffi::reaction_num_products(self.as_inner()))
    }

    /// Number of agent molecules (catalysts, solvents, …).
    pub fn num_agents(&self) -> u32 {
        crate::with_ob(|| ffi::reaction_num_agents(self.as_inner()))
    }

    /// The reactant at 0-based `index` as a standalone [`Molecule`] (a copy),
    /// or `None` if out of range.
    pub fn reactant(&self, index: u32) -> Option<Molecule> {
        let m = crate::with_ob(|| ffi::reaction_reactant(self.as_inner(), index));
        (!m.is_null()).then(|| Molecule::from_inner(m))
    }

    /// The product at 0-based `index` as a standalone [`Molecule`] (a copy),
    /// or `None` if out of range.
    pub fn product(&self, index: u32) -> Option<Molecule> {
        let m = crate::with_ob(|| ffi::reaction_product(self.as_inner(), index));
        (!m.is_null()).then(|| Molecule::from_inner(m))
    }

    /// The agent at 0-based `index` as a standalone [`Molecule`] (a copy), or
    /// `None` if out of range.
    pub fn agent(&self, index: u32) -> Option<Molecule> {
        let m = crate::with_ob(|| ffi::reaction_agent(self.as_inner(), index));
        (!m.is_null()).then(|| Molecule::from_inner(m))
    }

    /// Append a copy of `mol` to the reactant list.
    pub fn add_reactant(&mut self, mol: &Molecule) -> &mut Self {
        crate::with_ob(|| ffi::reaction_add_reactant(self.inner.pin_mut(), mol.as_inner()));
        self
    }

    /// Append a copy of `mol` to the product list.
    pub fn add_product(&mut self, mol: &Molecule) -> &mut Self {
        crate::with_ob(|| ffi::reaction_add_product(self.inner.pin_mut(), mol.as_inner()));
        self
    }

    /// Append a copy of `mol` to the agent list.
    pub fn add_agent(&mut self, mol: &Molecule) -> &mut Self {
        crate::with_ob(|| ffi::reaction_add_agent(self.inner.pin_mut(), mol.as_inner()));
        self
    }

    /// The reaction title (empty if unset).
    pub fn title(&self) -> String {
        crate::with_ob(|| ffi::reaction_title(self.as_inner()))
    }

    /// Set the reaction title.
    pub fn set_title(&mut self, title: &str) -> &mut Self {
        crate::with_ob(|| ffi::reaction_set_title(self.inner.pin_mut(), title));
        self
    }

    /// The reaction comment (empty if unset).
    pub fn comment(&self) -> String {
        crate::with_ob(|| ffi::reaction_comment(self.as_inner()))
    }

    /// Set the reaction comment.
    pub fn set_comment(&mut self, comment: &str) -> &mut Self {
        crate::with_ob(|| ffi::reaction_set_comment(self.inner.pin_mut(), comment));
        self
    }

    /// Whether the reaction is marked reversible.
    pub fn is_reversible(&self) -> bool {
        crate::with_ob(|| ffi::reaction_is_reversible(self.as_inner()))
    }

    /// Mark the reaction reversible (or not).
    pub fn set_reversible(&mut self, value: bool) -> &mut Self {
        crate::with_ob(|| ffi::reaction_set_reversible(self.inner.pin_mut(), value));
        self
    }

    fn as_inner(&self) -> &ffi::Reaction {
        self.inner.as_ref().expect("Reaction is never null")
    }
}

impl Default for Reaction {
    fn default() -> Self {
        Reaction::new()
    }
}
