#ifndef COOT_LIGAND_EDITOR_CANVAS_MODEL_HPP
#define COOT_LIGAND_EDITOR_CANVAS_MODEL_HPP
#include <cstddef>
#include <gtk/gtk.h>
#include <memory>
#include <vector>
#include <variant>
#include <optional>
#include <rdkit/GraphMol/RWMol.h>

namespace coot {
namespace ligand_editor_canvas {


/// For Edit Undo/Redo
class Operation {
    public:
    enum class OpType: unsigned char {
        ElementInsertion
    };
    class ElementInsertion {

    };
    private:
    union {
        ElementInsertion element_insertion;
    };
    OpType variant;
    public:
    
};

typedef std::vector<Operation> OperationStack;

/// Drawing-friendly representation of RDKit molecule
class CanvasMolecule {
    /// Used for lowering for RDKit.
    /// ABI must be kept compatible with BondType
    enum class BondTypeRaw: unsigned char {
        Single,
        Double,
        Triple,
        Aromatic
    };
    public:
    enum class AtomColor: unsigned char {
        /// Carbon and hydrogens
        Black,
        /// For Chlorine
        Green,
        /// For Nitrogen
        Blue,
        /// For Oxygen
        Red
        // are there more colors?
    };
    
    struct Atom {
        std::string symbol;
        AtomColor color;
        /// Position on canvas (x axis)
        float x;
        /// Position on canvas (y axis)
        float y;
        /// Corresponds to RDKit atom index
        unsigned int idx;
        bool highlighted;
    };
    enum class BondType: unsigned char {
        Single,
        Double,
        Triple
    };
    // todo: geometry support
    struct Bond {
        union {
            /// Used after lowering has been completed
            BondType type;
            /// Used while lowering from RDKit
            BondTypeRaw raw_type;
        };
        float first_atom_x;
        float first_atom_y;
        unsigned int first_atom_idx;
        float second_atom_x;
        float second_atom_y;
        unsigned int second_atom_idx;
        bool highlighted;
    };

    typedef std::optional<std::variant<CanvasMolecule::Atom,CanvasMolecule::Bond>> MaybeAtomOrBond;
    private:

    static const float ATOM_HITBOX_RADIUS;
    static const float BOND_DISTANCE_BOUNDARY;
    static const float BASE_SCALE_FACTOR;

    static BondTypeRaw bond_type_raw_from_rdkit(RDKit::Bond::BondType) noexcept;
    // static BondType bond_type_from_raw(BondTypeRaw);
    static AtomColor atom_color_from_rdkit(const RDKit::Atom *) noexcept;
    static std::tuple<float,float,float> atom_color_to_rgb(AtomColor) noexcept;
    static std::string atom_color_to_html(AtomColor) noexcept;

    std::shared_ptr<RDKit::RWMol> rdkit_molecule;
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;

    float _x_offset;
    float _y_offset;


    /// Computes the scale used for drawing
    /// And interfacing with screen coordinates
    float get_scale() const noexcept;

    public:

    CanvasMolecule(std::shared_ptr<RDKit::RWMol> rdkit_mol);

    /// Clears the drawing-friendly 2D representation data
    /// and re-creates it from the internal RDKit::RWMol
    void lower_from_rdkit();

    /// Changes the relative placement of the molecule on the screen
    void set_offset_from_bounds(const graphene_rect_t *bounds) noexcept;

    /// Draws the molecule on the widget.
    /// Be sure to call set_offset_from_bounds before drawing.
    void draw(GtkSnapshot* snapshot, PangoLayout* pango_layout, const graphene_rect_t *bounds) const noexcept;

    /// Checks if any object matches the click coordinates passed as arguments.
    /// Returns the thing that was clicked on (or nullopt if there's no match).
    MaybeAtomOrBond resolve_click(int x, int y) const noexcept;

    void highlight_atom(int atom_idx);
    void highlight_bond(int atom_a, int atom_b); 
    void clear_highlights();

};


} // namespace ligand_editor_canvas
} // namesapce coot

#endif // COOT_LIGAND_EDITOR_CANVAS_MODEL_HPP