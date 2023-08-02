#include "ligand_builder_ui.hpp"
#include "about_dialog.hpp"
#include "ligand_builder_state.hpp"
#include "ligand_editor_canvas.hpp"
#include "ligand_editor_canvas/core.hpp"
#include "ligand_editor_canvas/model.hpp"
#include "ligand_editor_canvas/tools.hpp"

void build_main_window(GtkWindow* win, CootLigandEditorCanvas* canvas, GtkLabel* status_label) {
    using namespace coot::ligand_editor_canvas;
    using BondModifierMode = coot::ligand_editor_canvas::BondModifier::BondModifierMode;
    using Element = coot::ligand_editor_canvas::ElementInsertion::Element;
    using Structure = coot::ligand_editor_canvas::StructureInsertion::Structure;
    using TransformMode = coot::ligand_editor_canvas::TransformManager::Mode;

    GtkWidget* mainbox = gtk_box_new(GTK_ORIENTATION_VERTICAL,5);

    gtk_window_set_child(win, mainbox);
    gtk_widget_set_margin_start(mainbox,10);
    gtk_widget_set_margin_end(mainbox,10);
    // Top toolbars

    GtkWidget* motions_toolbar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL,5);
    gtk_box_append(GTK_BOX(mainbox), motions_toolbar);
    GtkWidget* tools_toolbar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL,5);
    gtk_box_append(GTK_BOX(mainbox), tools_toolbar);
    GtkWidget* utils_toolbar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL,5);
    gtk_box_append(GTK_BOX(mainbox), utils_toolbar);

    GtkWidget* move_button = gtk_button_new_with_label("Move");
    g_signal_connect(move_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(TransformTool(TransformMode::Translation)));
    }), canvas);
    gtk_box_append(GTK_BOX(motions_toolbar), move_button);
    GtkWidget* rotate_button = gtk_button_new_with_label("Rotate");
    g_signal_connect(rotate_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(TransformTool(TransformMode::Rotation)));
    }), canvas);
    gtk_box_append(GTK_BOX(motions_toolbar), rotate_button);
    GtkWidget* flip_x_button = gtk_button_new_with_label("Flip around X");
    g_signal_connect(flip_x_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(FlipTool(FlipMode::Horizontal)));
    }), canvas);
    gtk_box_append(GTK_BOX(motions_toolbar), flip_x_button);
    GtkWidget* flip_y_button = gtk_button_new_with_label("Flip around Y");
    g_signal_connect(flip_y_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(FlipTool(FlipMode::Vertical)));
    }), canvas);
    gtk_box_append(GTK_BOX(motions_toolbar), flip_y_button);

    GtkWidget* single_bond_button = gtk_button_new_with_label("Single Bond");
    g_signal_connect(single_bond_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(BondModifier(BondModifierMode::Single)));
    }), canvas);
    gtk_box_append(GTK_BOX(tools_toolbar), single_bond_button);
    GtkWidget* double_bond_button = gtk_button_new_with_label("Double Bond");
    g_signal_connect(double_bond_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(BondModifier(BondModifierMode::Double)));
    }), canvas);
    gtk_box_append(GTK_BOX(tools_toolbar), double_bond_button);
    GtkWidget* triple_bond_button = gtk_button_new_with_label("Triple Bond");
    g_signal_connect(triple_bond_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(BondModifier(BondModifierMode::Triple)));
    }), canvas);
    gtk_box_append(GTK_BOX(tools_toolbar), triple_bond_button);
    GtkWidget* stereo_out_modifier_button = gtk_button_new_with_label("Geometry Tool");
    g_signal_connect(stereo_out_modifier_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(GeometryModifier()));
    }), canvas);
    gtk_box_append(GTK_BOX(tools_toolbar), stereo_out_modifier_button);
    GtkWidget* charge_modifier_button = gtk_button_new_with_label("Charge Tool");
    g_signal_connect(charge_modifier_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ChargeModifier()));
    }), canvas);
    gtk_box_append(GTK_BOX(tools_toolbar), charge_modifier_button);
    GtkWidget* delete_button = gtk_button_new_with_label("Delete");
    g_signal_connect(delete_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(DeleteTool()));
    }), canvas);
    gtk_box_append(GTK_BOX(tools_toolbar), delete_button);
    GtkWidget* format_button = gtk_button_new_with_label("Format");
    g_signal_connect(format_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(FormatTool()));
    }), canvas);
    gtk_box_append(GTK_BOX(tools_toolbar), format_button);
    GtkWidget* delete_hydrogens_button = gtk_button_new_with_label("Delete Hydrogens");
    g_signal_connect(delete_hydrogens_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(RemoveHydrogensTool()));
    }), canvas);
    gtk_box_append(GTK_BOX(tools_toolbar), delete_hydrogens_button);
    GtkWidget* smiles_button = gtk_button_new_with_label("SMILES");
    gtk_box_append(GTK_BOX(utils_toolbar), smiles_button);
    GtkWidget* buttom_EnvResidues = gtk_button_new_with_label("Env. Residues");
    gtk_box_append(GTK_BOX(utils_toolbar), buttom_EnvResidues);
    GtkWidget* buttom_Key = gtk_button_new_with_label("Key");
    gtk_box_append(GTK_BOX(utils_toolbar), buttom_Key);
    GtkWidget* info_button = gtk_button_new_with_label("Info");
    gtk_box_append(GTK_BOX(utils_toolbar), info_button);
    
    // Carbon ring picker
    GtkWidget* carbon_ring_picker = gtk_box_new(GTK_ORIENTATION_HORIZONTAL,5);
    gtk_box_append(GTK_BOX(mainbox), carbon_ring_picker);

    GtkWidget* button_3C = gtk_button_new_with_label("3-C");
    gtk_box_append(GTK_BOX(carbon_ring_picker), button_3C);
    g_signal_connect(button_3C, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(StructureInsertion(Structure::CycloPropaneRing)));

    }), canvas);
    GtkWidget* button_4C = gtk_button_new_with_label("4-C");
    gtk_box_append(GTK_BOX(carbon_ring_picker), button_4C);
    g_signal_connect(button_4C, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(StructureInsertion(Structure::CycloButaneRing)));

    }), canvas);
    GtkWidget* button_5C = gtk_button_new_with_label("5-C");
    gtk_box_append(GTK_BOX(carbon_ring_picker), button_5C);
    g_signal_connect(button_5C, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(StructureInsertion(Structure::CycloPentaneRing)));

    }), canvas);
    GtkWidget* button_6C = gtk_button_new_with_label("6-C");
    gtk_box_append(GTK_BOX(carbon_ring_picker), button_6C);
    g_signal_connect(button_6C, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(StructureInsertion(Structure::CycloHexaneRing)));

    }), canvas);
    GtkWidget* button_6Arom = gtk_button_new_with_label("6-Arom");
    gtk_box_append(GTK_BOX(carbon_ring_picker), button_6Arom);
    g_signal_connect(button_6Arom, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(StructureInsertion(Structure::BenzeneRing)));

    }), canvas);
    GtkWidget* button_7C = gtk_button_new_with_label("7-C");
    gtk_box_append(GTK_BOX(carbon_ring_picker), button_7C);
    g_signal_connect(button_7C, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(StructureInsertion(Structure::CycloHeptaneRing)));

    }), canvas);
    GtkWidget* button_8C = gtk_button_new_with_label("8-C");
    gtk_box_append(GTK_BOX(carbon_ring_picker), button_8C);
    g_signal_connect(button_8C, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(StructureInsertion(Structure::CycloOctaneRing)));

    }), canvas);

    GtkWidget* canvas_space = gtk_box_new(GTK_ORIENTATION_HORIZONTAL,0);
    gtk_box_append(GTK_BOX(mainbox), canvas_space);
    // Canvas space: chemical element picker
    GtkWidget* chem_element_picker = gtk_box_new(GTK_ORIENTATION_VERTICAL,5);
    gtk_widget_set_margin_top(chem_element_picker,10);

    gtk_box_append(GTK_BOX(canvas_space), chem_element_picker);

    GtkWidget* C_button = gtk_button_new_with_label("C");
    g_signal_connect(C_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::C)));
    }), canvas);
    gtk_box_append(GTK_BOX(chem_element_picker), C_button);
    GtkWidget* N_button = gtk_button_new_with_label("N");
    g_signal_connect(N_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::N)));
    }), canvas);
    gtk_box_append(GTK_BOX(chem_element_picker), N_button);
    GtkWidget* O_button = gtk_button_new_with_label("O");
    g_signal_connect(O_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::O)));
    }), canvas);
    gtk_box_append(GTK_BOX(chem_element_picker), O_button);
    GtkWidget* S_button = gtk_button_new_with_label("S");
    g_signal_connect(S_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::S)));
    }), canvas);
    gtk_box_append(GTK_BOX(chem_element_picker), S_button);
    GtkWidget* P_button = gtk_button_new_with_label("P");
    g_signal_connect(P_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::P)));
    }), canvas);
    gtk_box_append(GTK_BOX(chem_element_picker), P_button);
    GtkWidget* H_button = gtk_button_new_with_label("H");
    g_signal_connect(H_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::H)));
    }), canvas);
    gtk_box_append(GTK_BOX(chem_element_picker), H_button);
    GtkWidget* F_button = gtk_button_new_with_label("F");
    g_signal_connect(F_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::F)));
    }), canvas);
    gtk_box_append(GTK_BOX(chem_element_picker), F_button);
    GtkWidget* Cl_button = gtk_button_new_with_label("Cl");
    g_signal_connect(Cl_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::Cl)));
    }), canvas);
    gtk_box_append(GTK_BOX(chem_element_picker), Cl_button);
    GtkWidget* Br_button = gtk_button_new_with_label("Br");
    g_signal_connect(Br_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::Br)));
    }), canvas);
    gtk_box_append(GTK_BOX(chem_element_picker), Br_button);
    GtkWidget* I_button = gtk_button_new_with_label("I");
    g_signal_connect(I_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::I)));
    }), canvas);
    gtk_box_append(GTK_BOX(chem_element_picker), I_button);
    GtkWidget* X_button = gtk_button_new_with_label("X");
    g_signal_connect(X_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        ((coot::ligand_editor::LigandBuilderState*)user_data)->run_choose_element_dialog();
    }), nullptr);
}

void coot::ligand_editor::setup_actions(coot::ligand_editor::LigandBuilderState* state, GtkApplicationWindow* win, GtkBuilder* builder) {
    auto new_action = [win](const char* action_name, GCallback func, gpointer userdata = nullptr){
        std::string detailed_action_name = "win.";
        detailed_action_name += action_name;
        GSimpleAction* action = g_simple_action_new(action_name,nullptr);
        g_action_map_add_action(G_ACTION_MAP(win), G_ACTION(action));
        g_signal_connect(action, "activate", func, userdata);
        //return std::make_pair(detailed_action_name,action);
    };

    auto new_stateful_action = [win](const char* action_name,const GVariantType *state_type, GVariant* default_state, GCallback func, gpointer userdata = nullptr){
        std::string detailed_action_name = "win.";
        detailed_action_name += action_name;
        GSimpleAction* action = g_simple_action_new_stateful(action_name, state_type, default_state);
        g_action_map_add_action(G_ACTION_MAP(win), G_ACTION(action));
        g_signal_connect(action, "activate", func, userdata);
        //return std::make_pair(detailed_action_name,action);
    };

    using ExportMode = coot::ligand_editor::LigandBuilderState::ExportMode;

    // File
    new_action("file_new", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LigandBuilderState*)user_data)->file_new();
    }),state);
    new_action("file_open", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LigandBuilderState*)user_data)->file_open();
    }),state);
    new_action("import_from_smiles", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LigandBuilderState*)user_data)->load_from_smiles();
    }),state);
    new_action("import_molecule", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LigandBuilderState*)user_data)->file_import_molecule();
    }),state);
    new_action("fetch_molecule", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LigandBuilderState*)user_data)->file_fetch_molecule();
    }),state);
    new_action("file_save", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LigandBuilderState*)user_data)->file_save();
    }),state);
    new_action("file_save_as", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LigandBuilderState*)user_data)->file_save_as();
    }),state);
    new_action("export_pdf", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LigandBuilderState*)user_data)->file_export(ExportMode::PDF);
    }),state);
    new_action("export_png", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LigandBuilderState*)user_data)->file_export(ExportMode::PNG);
    }),state);
    new_action("export_svg", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LigandBuilderState*)user_data)->file_export(ExportMode::SVG);
    }),state);
    new_action("file_exit", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LigandBuilderState*)user_data)->file_exit();
    }),state);

    // Edit;
    new_action("undo", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LigandBuilderState*)user_data)->edit_undo();
    }),state);
    new_action("redo", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LigandBuilderState*)user_data)->edit_redo();
    }),state);
    // Display

    using coot::ligand_editor_canvas::DisplayMode;
    GVariant* display_mode_action_defstate = g_variant_new("s",coot::ligand_editor_canvas::display_mode_to_string(DisplayMode::Standard));
    new_stateful_action(
        "switch_display_mode", 
        G_VARIANT_TYPE_STRING,
        display_mode_action_defstate, 
        G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
            const gchar* mode_name = g_variant_get_string(parameter,nullptr);
            auto mode = coot::ligand_editor_canvas::display_mode_from_string(mode_name);
            if(mode.has_value()) {
                ((LigandBuilderState*)user_data)->switch_display_mode(mode.value());
                g_simple_action_set_state(self, parameter);
            } else {
                g_error("Could not parse display mode from string!: '%s'",mode_name);
            }
        }
    ),state);

    // Help

    new_action("show_about_dialog", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        auto* about_dialog = GTK_WINDOW(user_data);
        gtk_window_present(GTK_WINDOW(about_dialog));
    }),gtk_builder_get_object(builder, "layla_about_dialog"));

}