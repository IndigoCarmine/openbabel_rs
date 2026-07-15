//! SMARTS-based structural transformations (SMIRKS-like).

use cxx::UniquePtr;
use openbabel_sys::ffi;

use crate::error::Error;
use crate::mol::Molecule;

/// A compiled transformation that rewrites substructures matching a reactant
/// SMARTS into a product SMARTS — OpenBabel's `OBChemTsfm`, the mechanism
/// behind its pH model. Atoms are carried from reactant to product by their
/// SMARTS class labels (`:n`); the transform can change elements, charges,
/// bond orders, and add or delete atoms/bonds.
///
/// ```no_run
/// use openbabel::{Molecule, Transform};
/// // Deprotonate a carboxylic acid: -C(=O)OH -> -C(=O)[O-]
/// let t = Transform::new("[C:1](=O)[OH:2]", "[C:1](=O)[O-:2]").unwrap();
/// let mut mol = Molecule::parse("CC(=O)O", "smi").unwrap(); // acetic acid
/// assert!(t.apply(&mut mol));
/// assert_eq!(mol.total_charge(), -1);
/// ```
pub struct Transform {
    inner: UniquePtr<ffi::Transform>,
}

impl Transform {
    /// Compile a transformation from a reactant SMARTS to a product SMARTS.
    ///
    /// Returns [`Error::InvalidTransform`] if either pattern is malformed.
    pub fn new(reactant_smarts: &str, product_smarts: &str) -> Result<Self, Error> {
        let inner = crate::with_ob(|| ffi::transform_new(reactant_smarts, product_smarts));
        if inner.is_null() {
            Err(Error::InvalidTransform {
                reactant: reactant_smarts.to_string(),
                product: product_smarts.to_string(),
            })
        } else {
            Ok(Transform { inner })
        }
    }

    /// Apply the transformation to every match in `mol`, editing it in place.
    /// Returns `false` if nothing matched (the molecule is then unchanged).
    pub fn apply(&self, mol: &mut Molecule) -> bool {
        crate::with_ob(|| ffi::transform_apply(self.as_inner(), mol.as_inner_pin_mut()))
    }

    fn as_inner(&self) -> &ffi::Transform {
        self.inner.as_ref().expect("Transform is never null")
    }
}
