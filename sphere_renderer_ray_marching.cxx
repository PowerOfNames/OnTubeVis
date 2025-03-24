#include "sphere_renderer_ray_marching.h"


#include <cgv_gl/gl/gl.h>
#include <cgv_gl/gl/gl_tools.h>

namespace cgv {
	namespace render {
		sphere_renderer_ray_marching& ref_sphere_renderer_ray_marching(context& ctx, int ref_count_change)
		{
			static int ref_count = 0;
			static sphere_renderer_ray_marching r;
			r.manage_singleton(ctx, "sphere_renderer_ray_marching", ref_count, ref_count_change);
			return r;
		}

		render_style* sphere_renderer_ray_marching::create_render_style() const
		{
			return new sphere_render_ray_marching_style();
		}

		sphere_render_ray_marching_style::sphere_render_ray_marching_style()
		{
			radius_scale = 1;
			radius = 1;
			use_group_radius = false;
			morph_sin_factor = 0.0f;

			blend_width_in_pixel = 0.0f;

			//Ray marching defines:
			rm.epsilon = 0.003f;
			rm.max_iterations = 30;
			rm.fdg_delta = 0.05f;
			rm.show_bounding = false;

			glyph_color_mapping = vec4(0.0f, 0.1f, 0.0f, 0.1f);
		}

		sphere_renderer_ray_marching::sphere_renderer_ray_marching()
		{
			has_radii = false;
			has_group_radii = false;
			cull_per_primitive = false;
		}

		/// call this before setting attribute arrays to manage attribute array in given manager
		void sphere_renderer_ray_marching::enable_attribute_array_manager(const context& ctx, attribute_array_manager& aam)
		{
			surface_renderer::enable_attribute_array_manager(ctx, aam);
			if (has_attribute(ctx, "radius"))
				has_radii = true;
		}
		/// call this after last render/draw call to ensure that no other users of renderer change attribute arrays of given manager
		void sphere_renderer_ray_marching::disable_attribute_array_manager(const context& ctx, attribute_array_manager& aam)
		{
			surface_renderer::disable_attribute_array_manager(ctx, aam);
			has_radii = false;
		}
		void sphere_renderer_ray_marching::remove_radius_array(const context& ctx) {
			has_radii = false;
			remove_attribute_array(ctx, "radius");
		}
		///
		void sphere_renderer_ray_marching::set_y_view_angle(float _y_view_angle)
		{
			y_view_angle = _y_view_angle;
		}

		void sphere_renderer_ray_marching::update_defines(shader_define_map& defines)
		{
			const sphere_render_ray_marching_style& rs = get_style<sphere_render_ray_marching_style>();
			defines.clear();

			shader_code::set_define(defines, "RM_EPSILON", rs.rm.epsilon, 0.001f);
			shader_code::set_define(defines, "RM_MAX_ITERATIONS", rs.rm.max_iterations, (uint32_t)30);
			shader_code::set_define(defines, "RM_FDG_DELTA", rs.rm.fdg_delta, 0.0005f);
			shader_code::set_define(defines, "RM_SHOW_BOUNDING", rs.rm.show_bounding, false);
		}

		bool sphere_renderer_ray_marching::build_shader_program(context& ctx, shader_program& prog, const shader_define_map& defines)
		{
			return prog.build_program(ctx, "sphere_ray_marched.glpr", true, defines);
		}
		bool sphere_renderer_ray_marching::validate_attributes(const context& ctx) const
		{
			const sphere_render_ray_marching_style& srs = get_style<sphere_render_ray_marching_style>();
			bool res = surface_renderer::validate_attributes(ctx);
			if (!has_group_radii && srs.use_group_radius) {
				ctx.error("sphere_renderer_ray_marching::validate_attributes() group_radii not set");
				res = false;
			}
			return res;
		}
		bool sphere_renderer_ray_marching::enable(context& ctx)
		{
			const sphere_render_ray_marching_style& rs = get_style<sphere_render_ray_marching_style>();

			if (!surface_renderer::enable(ctx))
				return false;

			if (!ref_prog().is_linked())
				return false;

			if (!has_radii)
				ref_prog().set_attribute(ctx, "radius", rs.radius);

			ref_prog().set_uniform(ctx, "use_group_radius", rs.use_group_radius);
			ref_prog().set_uniform(ctx, "radius_scale", rs.radius_scale);

			ref_prog().set_uniform(ctx, "glyph_cc_param", rs.glyph_color_mapping);

			float pixel_extent_per_depth = (float)(2.0 * tan(0.5 * 0.0174532925199 * y_view_angle) / ctx.get_height());
			ref_prog().set_uniform(ctx, "pixel_extent_per_depth", pixel_extent_per_depth);
			ref_prog().set_uniform(ctx, "blend_width_in_pixel", rs.blend_width_in_pixel);
			ref_prog().set_uniform(ctx, "morph_sin_factor", rs.morph_sin_factor);

			return true;
		}

		bool sphere_renderer_ray_marching::disable(context& ctx)
		{
			const sphere_render_ray_marching_style& srs = get_style<sphere_render_ray_marching_style>();

			if (!attributes_persist()) {
				has_radii = false;
				has_group_radii = false;
			}
			return surface_renderer::disable(ctx);
		}

		void sphere_renderer_ray_marching::draw(context& ctx, size_t start, size_t count, bool use_strips, bool use_adjacency, uint32_t strip_restart_index)
		{
			draw_impl(ctx, PT_POINTS, start, count, false, false, -1);
		}

		bool sphere_render_ray_marching_style_reflect::self_reflect(cgv::reflect::reflection_handler& rh)
		{
			return
				rh.reflect_base(*static_cast<surface_render_style*>(this)) &&
				rh.reflect_member("radius", radius) &&
				rh.reflect_member("use_group_radius", use_group_radius) &&
				rh.reflect_member("radius_scale", radius_scale) &&
				rh.reflect_member("blend_width_in_pixel", blend_width_in_pixel) &&
				rh.reflect_member("morph_sin_factor", morph_sin_factor) &&
				rh.reflect_member("rm.epsilon", rm.epsilon) &&
				rh.reflect_member("rm.max_iterations", rm.max_iterations) &&
				rh.reflect_member("rm.show_bounding", rm.show_bounding) &&
				rh.reflect_member("rm.fdg_delta", rm.fdg_delta);
		}

		cgv::reflect::extern_reflection_traits<sphere_render_ray_marching_style, sphere_render_ray_marching_style_reflect> get_reflection_traits(const sphere_render_ray_marching_style&)
		{
			return cgv::reflect::extern_reflection_traits<sphere_render_ray_marching_style, sphere_render_ray_marching_style_reflect>();
		}
	}
}
#include <cgv/gui/provider.h>


namespace cgv {
	namespace gui {

		struct sphere_render_ray_marching_style_gui_creator : public gui_creator {
			/// attempt to create a gui and return whether this was successful
			bool create(provider* p, const std::string& label,
				void* value_ptr, const std::string& value_type,
				const std::string& gui_type, const std::string& options, bool*) {
				if (value_type != cgv::type::info::type_name<cgv::render::sphere_render_ray_marching_style>::get_name())
					return false;
				cgv::render::sphere_render_ray_marching_style* rs_ptr = reinterpret_cast<cgv::render::sphere_render_ray_marching_style*>(value_ptr);
				cgv::base::base* b = dynamic_cast<cgv::base::base*>(p);

				if (p->begin_tree_node("Sphere Ray Marching Parameters", rs_ptr->rm, false)) {
					p->align("\a");
					p->add_member_control(b, "Epsilon", rs_ptr->rm.epsilon, "value_slider", "min=0;max=0.1;step=0.001;log=true");
					p->add_member_control(b, "Max Iterations", rs_ptr->rm.max_iterations, "value_slider", "min=0;max=50;step=1;unsigned=true");
					p->add_member_control(b, "FinDiffGrad Delta", rs_ptr->rm.fdg_delta, "value_slider", "min=0.01;max=1.0;step=0.01;ticks=true");
					p->add_member_control(b, "FinDiffGrad Delta", rs_ptr->rm.fdg_delta, "value_slider", "min=0.01;max=1.0;step=0.01;ticks=true");
					p->add_member_control(b, "Show Bounding", rs_ptr->rm.show_bounding, "check");
					p->align("\b");
					p->end_tree_node(rs_ptr->rm);
				}
				p->add_member_control(b, "Morph Sin Factor", rs_ptr->morph_sin_factor, "value_slider", "min=0;max=100.0;step=0.1;ticks=true");

				return true;
			}
		};

#include <cgv_gl/gl/lib_begin.h>

		cgv::gui::gui_creator_registration<sphere_render_ray_marching_style_gui_creator> sphere_render_ray_marching_rs_gc_reg("sphere_render_ray_marching_style_gui_creator");
	}
}
