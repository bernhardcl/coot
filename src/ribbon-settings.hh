

namespace coot {

   // ------------ ribbons settings (should be rather general?ccp4mg? FIXME) ------------
   //

   class ribbon_settings_t {
   public:
      int solid_quality; // fast, smooth, delux (0,1,2) can we use?
      float cylinder_width; // needed?
      float ribbon_width;
      float alpha_helix_width;
      float arrow_width;
      float arrow_length;
      float worm_width;
      float helix_tube_diameter;
      int ribbon_style; // oval, flat, round (0,1,2)
      int two_colour_ribbon; // normal, inside grey (0,1)
      int helix_style;  // oval, flat, flat/round, fancy (0,1,2,3)
      int flatten_loop; //shorten loops
      int flatten_beta; // smooth b-sheets
      int spline_beta_flat; // needed ?
      ribbon_settings_t() {
         solid_quality = 1;
         cylinder_width = 0.2;
         ribbon_width = 2.5;
         alpha_helix_width = 1.5;
         arrow_width = 2.2;
         arrow_length = 2.2;
         worm_width = 0.2;
         helix_tube_diameter = 1.5;
         ribbon_style = 0;
         two_colour_ribbon = 0;
         helix_style = 0;
         flatten_loop = 0;
         flatten_beta = 1;
         spline_beta_flat = 1;
      }
   };



}
