#pragma once

#include <cgv_gl/surface_renderer.h>

#include <cgv_gl/gl/lib_begin.h>

namespace cgv { // @<
	namespace render { // @<
		class cone_renderer_ray_marching;

		//! reference to a singleton cone renderer that is shared among drawables
		/*! the second parameter is used for reference counting. Use +1 in your init method,
			-1 in your clear method and default 0 argument otherwise. If internal reference
			counter decreases to 0, singleton renderer is destructed. */

		extern  cone_renderer_ray_marching& ref_cone_renderer_ray_marching(context& ctx, int ref_count_change = 0);
		struct cone_render_ray_marching_style : public surface_render_style
		{
			/// multiplied to the sphere radii, initialized to 1
			float radius_scale;

			/// default value assigned to radius attribute in \c enable method of cone renderer, set to 1 in constructor
			float radius;

			bool show_caps;
			bool rounded_caps;
			bool composite_arrow;

			vec4 glyph_color_mapping;

			struct RayMarching
			{
				float epsilon;
				uint32_t max_iterations;
				float fdg_delta;
				bool show_bounding;

			}rm;

			/// construct with default values
			cone_render_ray_marching_style();
		};

		/// renderer that supports raycasting of cones
		class cone_renderer_ray_marching : public surface_renderer
		{
		protected:
			bool has_radii;
			/// the shader defines used to build the shader, used to comapre against new defines to determine if the shader needs to be rebuilt
			shader_define_map shader_defines;
			/// overload to allow instantiation of cone_renderer_ray_marching
			render_style* create_render_style() const;
			/// update shader defines based on render style
			void update_defines(shader_define_map& defines);
			/// build cone program
			bool build_shader_program(context& ctx, shader_program& prog, const shader_define_map& defines);

		public:
			/// initializes member variables
			cone_renderer_ray_marching();
			/// call this before setting attribute arrays to manage attribute array in given manager
			void enable_attribute_array_manager(const context& ctx, attribute_array_manager& aam);
			/// call this after last render/draw call to ensure that no other users of renderer change attribute arrays of given manager
			void disable_attribute_array_manager(const context& ctx, attribute_array_manager& aam);
			///
			bool enable(context& ctx);
			///
			template <typename T = float>
			void set_radius_array(const context& ctx, const std::vector<T>& radii) { has_radii = true; set_attribute_array(ctx, "radius", radii); }
			///
			template <typename T = float>
			void set_radius_array(const context& ctx, const T* radii, size_t nr_elements, unsigned stride_in_bytes = 0) { has_radii = true; set_attribute_array(ctx, "radius", radii, nr_elements, stride_in_bytes); }
			/// remove the radius attribute
			void remove_radius_array(const context& ctx);

			/// use this function if you store spheres in vec4 with the 4th component the radius
			template <typename T = float>
			void set_sphere_array(const context& ctx, const std::vector<cgv::math::fvec<T, 4> >& spheres) {
				set_composed_attribute_array(ctx, "position", &spheres.front(), spheres.size(), reinterpret_cast<const cgv::math::fvec<T, 3>&>(spheres.front()));
				ref_composed_attribute_array(ctx, "radius", "position", &spheres.front(), spheres.size(), spheres[0][3]);
				has_positions = true;
				has_radii = true;
			}
			///
			bool validate_attributes(const context& ctx) const;
			///
			bool disable(context& ctx);
			/// convenience function to render with default settings
			void draw(context& ctx, size_t start, size_t count,
				bool use_strips = false, bool use_adjacency = false, uint32_t strip_restart_index = -1);
			/// the clear function destructs the shader program and resets the texture pointers
			virtual void clear(const context& ctx);
		};

		struct cone_render_ray_marching_style_reflect : public cone_render_ray_marching_style
		{
			bool self_reflect(cgv::reflect::reflection_handler& rh);
		};

		extern cgv::reflect::extern_reflection_traits<cone_render_ray_marching_style, cone_render_ray_marching_style_reflect> get_reflection_traits(const cone_render_ray_marching_style&);
	}
}

#include <cgv/config/lib_end.h>