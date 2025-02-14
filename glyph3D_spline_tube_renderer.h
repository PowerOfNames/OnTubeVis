#pragma once
#include "cgv_gl/surface_renderer.h"
#include <cgv_reflect_types/media/color.h>


#include "cgv_gl/gl/lib_begin.h"
#include "textured_spline_tube_renderer.h"

namespace cgv { // @<
	namespace render { // @<
		class glyph3D_spline_tube_renderer;


		extern glyph3D_spline_tube_renderer& ref_glyph3D_spline_tube_renderer(context& ctx, int ref_count_change = 0, textured_spline_tube_render_style* textured_rs_ptr = nullptr);

		/*!	Style to control the look of glyph3D spline tubes. */
		struct glyph3D_spline_tube_render_style : public surface_render_style
		{
			enum GlyphMethod
			{
				GM_TUBE_SURFACE_PIPELINE = 0,
				GM_BILLBOARD_PIPELINE = 1
			} glyph_method;

			/// construct with default values
			glyph3D_spline_tube_render_style();

			/// check wether chosen glyph medhod is textured
			inline bool is_tube_surface_pipeline(void) const {
				return glyph_method == 0;
			}
			/// check wether glyph method is billboard
			inline bool is_billboard_pipeline(void) const {
				return glyph_method == 1;
			}

			struct tube_backside_render_parameters {
				bool use_distance_check;
				bool clip_caps;
				float distance_tolerance_factor;
			};
			tube_backside_render_parameters tube_backside_render;

			bool cull_frontface;
			float front_transparency;
			float fresnel_reflection;
		};

		/// renderer that supports textured cubic hermite spline tubes
		//class glyph3D_spline_tube_renderer : public surface_renderer
		class glyph3D_spline_tube_renderer : public surface_renderer
		{
		protected:
			/// whether node ids are specified
			bool has_node_ids;
			/// whether radii are specified
			bool has_radii;
			/// whether tangents are specified
			bool has_tangents;
			/// position of the cyclopian eye point (differs from eye point in case of stereoscopic rendering)
			vec3 cyclopic_eye;
			/// camera view direction
			vec3 view_dir;
			/// viewport rectangle (offset and size)
			vec4 viewport;
			/// additional defines not dependant on the style and set from outside the renderer
			shader_define_map additional_defines;
			/// keep track of which line primitive was active the last time the renderer drew something


			//THESIS2:
			glyph3D_spline_tube_render_style::GlyphMethod last_active_glyph_method;
			//pointer to base textured_spline_tube_renderer_pointer to get access to standard tube rendering fields without copy and updating everything twice.
			textured_spline_tube_render_style* textured_rs;

			/// overload to allow instantiation of box_renderer
			render_style* create_render_style() const;
			/// update shader defines based on render style
			void update_defines(shader_define_map& defines);
			/// build rounded cone program
			bool build_shader_program(context& ctx, shader_program& prog, const shader_define_map& defines);

		public:
			/// initializes position_is_center to true 
			glyph3D_spline_tube_renderer();

			//THESIS2:
			inline void set_textured_spline_tube_render_style_ptr(textured_spline_tube_render_style* textured_rs_ptr) { textured_rs = textured_rs_ptr; }

			/// call this before setting attribute arrays to manage attribute array in given manager
			void enable_attribute_array_manager(const context& ctx, attribute_array_manager& aam);
			/// call this after last render/draw call to ensure that no other users of renderer change attribute arrays of given manager
			void disable_attribute_array_manager(const context& ctx, attribute_array_manager& aam);
			///
			void set_cyclopic_eye(const vec3& cyclopic_eye_pos) { this->cyclopic_eye = cyclopic_eye_pos; }
			///
			void set_view_dir(const vec3& view_dir) { this->view_dir = view_dir; }
			///
			void set_viewport(const vec4& viewport) { this->viewport = viewport; }
			/// set additional defines that do not depend on the style
			void set_additional_defines(shader_define_map& defines);
			///
			template <typename T = float>
			void set_node_id_array(const context& ctx, const std::vector<T>& node_ids) { has_node_ids = true; set_attribute_array(ctx, "node_ids", node_ids); }
			/// 
			template <typename T = float>
			void set_node_id_array(const context& ctx, const T* node_ids, size_t nr_elements, unsigned stride_in_bytes = 0) { has_node_ids = true; set_attribute_array(ctx, "node_ids", node_ids, nr_elements, stride_in_bytes); }
			///
			template <typename T = float>
			void set_radius_array(const context& ctx, const std::vector<T>& radii) { has_radii = true; set_attribute_array(ctx, "radius", radii); }
			/// 
			template <typename T = float>
			void set_radius_array(const context& ctx, const T* radii, size_t nr_elements, unsigned stride_in_bytes = 0) { has_radii = true; set_attribute_array(ctx, "radius", radii, nr_elements, stride_in_bytes); }
			///
			template <typename T = float>
			void set_tangent_array(const context& ctx, const std::vector<T>& tangents) { has_tangents = true; set_attribute_array(ctx, "tangent", tangents); }
			/// 
			template <typename T = float>
			void set_tangent_array(const context& ctx, const T* tangents, size_t nr_elements, unsigned stride_in_bytes = 0) { has_tangents = true; set_attribute_array(ctx, "tangent", tangents, nr_elements, stride_in_bytes); }
			///
			bool validate_attributes(const context& ctx) const;
			/// 
			bool enable(context& ctx);
			///
			bool disable(context& ctx);
			///
			void draw(context& ctx, size_t sunt, size_t count,
				bool use_strips = false, bool use_adjacency = false, uint32_t strip_resunt_index = -1);
		};

		struct glyph3D_spline_tube_render_style_reflect : public glyph3D_spline_tube_render_style
		{
			bool self_reflect(cgv::reflect::reflection_handler& rh);
		};
		extern cgv::reflect::extern_reflection_traits<glyph3D_spline_tube_render_style, glyph3D_spline_tube_render_style_reflect> get_reflection_traits(const glyph3D_spline_tube_render_style&);
	}
}

namespace cgv {
	namespace reflect {
		enum_reflection_traits<cgv::render::glyph3D_spline_tube_render_style::GlyphMethod> get_reflection_traits(const cgv::render::glyph3D_spline_tube_render_style::GlyphMethod&);
	}
}

#include <cgv/config/lib_end.h>


