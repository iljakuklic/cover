//! Maintain mapping between structured values and unique monotonically
//! increasing identifiers.

use std::collections::HashMap;

/// Automatically allocate numerical IDs to objects.
///
/// Upon insertion of an intem into this data structure, it gets assigned
/// a uniwue ID. Duplicate insertions give the same ID. IDs are assigned
/// in a monotonically increasing fashion, starting at 0.
struct IDMap<ObjT, ID> {
    fwd: HashMap<ObjT, ID>,
    bwd: Vec<ObjT>,
}

impl<ObjT, ID> IDMap<ObjT, ID>
where ObjT: Clone + Eq + std::hash::Hash, ID: Copy + From<usize> + Into<usize> {

    /// Create a new empty IDMap.
    pub fn new() -> IDMap<ObjT, ID> {
        IDMap { fwd: HashMap::new(), bwd: Vec::new() }
    }

    /// Lookup item but do not assign ID if the item is not already present.
    pub fn lookup(&self, obj: &ObjT) -> Option<ID> {
        self.fwd.get(obj).copied()
    }

    /// Get item ID, allocating a new ID if required.
    pub fn get(&mut self, obj: &ObjT) -> ID {
        if let Some(id) = self.fwd.get(&obj) {
            return *id;
        }
        let id: ID = self.bwd.len().into();
        self.fwd.insert(obj.clone(), id);
        self.bwd.push(obj.clone());
        id
    }

    /// Lookup item by ID. Panics if ID is not already present.
    pub fn back(&self, id: ID) -> &ObjT {
        &self.bwd[id.into()]
    }

    /// Iterate over all the IDs and corresponding items.
    pub fn items(&self) -> impl Iterator<Item=(ID, &ObjT)> {
        self.bwd.iter().enumerate().map(|(id, o)| (id.into(), o))
    }

    /// Length, i.e. number of unique items inserted.
    pub fn len(&self) -> usize {
        self.bwd.len()
    }
}

