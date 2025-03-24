#include "ellipsoid_renderer_ray_marching.h"


#include <cgv_gl/gl/gl.h>
#include <cgv_gl/gl/gl_tools.h>

namespace cgv {
	namespace render {
		ellipsoid_renderer_ray_marching& ref_ellipsoid_renderer_ray_marching(context& ctx, int ref_count_change)
		{
			static int ref_count = 0;
			static ellipsoid_renderer_ray_marching r;
			r.manage_singleton(ctx, "ellipsoid_renderer_ray_marching", ref_count, ref_count_change);
			return r;
		}

		render_style* ellipsoid_renderer_ray_marching::create_render_style() const
		{
			return new ellipsoid_render_ray_marching_style();
		}

		ellipsoid_render_ray_marching_style::ellipsoid_render_ray_marching_style()
		{
			size_scale = 1;
			size = 1;


			//Ray marching defines:
			rm.epsilon = 0.003f;
			rm.max_iterations = 30;
			rm.fdg_delta = 0.05f;
			rm.show_bounding = false;
			glyph_color_mapping = vec4(0.0f, 0.1f, 0.0f, 0.1f);
			
		}

		ellipsoid_renderer_ray_marching::ellipsoid_renderer_ray_marching()
		{
			has_sizes = false;
			has_orientations = false;
			cull_per_primitive = false;
		}
		/// call this before setting attribute arrays to manage attribute array in given manager
		void ellipsoid_renderer_ray_marching::enable_attribute_array_manager(const context& ctx, attribute_array_manager& aam)
		{
			surface_renderer::enable_attribute_array_manager(ctx, aam);
			if (has_attribute(ctx, "size"))
				has_sizes = true;
			if (has_attribute(ctx, "orientation"))
				has_orientations = true;
		}
		/// call this after last render/draw call to ensure that no other users of renderer change attribute arrays of given manager
		void ellipsoid_renderer_ray_marching::disable_attribute_array_manager(const context& ctx, attribute_array_manager& aam)
		{
			surface_renderer::disable_attribute_array_manager(ctx, aam);
			has_sizes = false;
			has_orientations = false;
		}
		void ellipsoid_renderer_ray_marching::remove_size_array(const context& ctx) {
			has_sizes = false;
			remove_attribute_array(ctx, "size");
		}

		void ellipsoid_renderer_ray_marching::update_defines(shader_define_map& defines)
		{
			const ellipsoid_render_ray_marching_style& rs = get_style<ellipsoid_render_ray_marching_style>();
			defines.clear();

			shader_code::set_define(defines, "RM_EPSILON", rs.rm.epsilon, 0.001f);
			shader_code::set_define(defines, "RM_MAX_ITERATIONS", rs.rm.max_iterations, (uint32_t)30);
			shader_code::set_define(defines, "RM_FDG_DELTA", rs.rm.fdg_delta, 0.0005f);
			shader_code::set_define(defines, "RM_SHOW_BOUNDING", rs.rm.show_bounding, false);
		}

		void ellipsoid_renderer_ray_marching::remove_orientation_array(const context& ctx) {
			has_orientations = false;
			remove_attribute_array(ctx, "orientation");
		}
		bool ellipsoid_renderer_ray_marching::build_shader_program(context& ctx, shader_program& prog, const shader_define_map& defines)
		{
			return prog.build_program(ctx, "ellipsoid_ray_marched.glpr", true, defines);
		}
		bool ellipsoid_renderer_ray_marching::validate_attributes(const context& ctx) const
		{
			const ellipsoid_render_ray_marching_style& rs = get_style<ellipsoid_render_ray_marching_style>();
			return surface_renderer::validate_attributes(ctx);
		}
		bool ellipsoid_renderer_ray_marching::enable(context& ctx)
		{
			const ellipsoid_render_ray_marching_style& rs = get_style<ellipsoid_render_ray_marching_style>();

			if (!surface_renderer::enable(ctx))
				return false;

			if (!ref_prog().is_linked())
				return false;

			if (!has_sizes)
				ref_prog().set_attribute(ctx, "size", rs.size);

			if (!has_orientations)
				ref_prog().set_attribute(ctx, "orientation", cgv::math::quaternion<float>());

			ref_prog().set_uniform(ctx, "size_scale", rs.size_scale);			
			ref_prog().set_uniform(ctx, "glyph_cc_param", rs.glyph_color_mapping);

			return true;
		}

		bool ellipsoid_renderer_ray_marching::disable(context& ctx)
		{
			if (!attributes_persist()) {
				has_sizes = false;
				has_orientations = false;
			}
			return surface_renderer::disable(ctx);
		}

		void ellipsoid_renderer_ray_marching::draw(context& ctx, size_t start, size_t count, bool use_strips, bool use_adjacency, uint32_t strip_restart_index)
		{
			draw_impl(ctx, PT_POINTS, start, count, false, false, -1);
		}

		bool ellipsoid_render_ray_marching_style_reflect::self_reflect(cgv::reflect::reflection_handler& rh)
		{
			return
				rh.reflect_base(*static_cast<surface_render_style*>(this)) &&
				rh.reflect_member("size_scale", size_scale) &&
				rh.reflect_member("rm.epsilon", rm.epsilon) &&
				rh.reflect_member("rm.max_iterations", rm.max_iterations) &&
				rh.reflect_member("rm.show_bounding", rm.show_bounding) &&
				rh.reflect_member("rm.fdg_delta", rm.fdg_delta);
		}

		cgv::reflect::extern_reflection_traits<ellipsoid_render_ray_marching_style, ellipsoid_render_ray_marching_style_reflect> get_reflection_traits(const ellipsoid_render_ray_marching_style&)
		{
			return cgv::reflect::extern_reflection_traits<ellipsoid_render_ray_marching_style, ellipsoid_render_ray_marching_style_reflect>();
		}
	}
}

#include <cgv/gui/provider.h>

namespace cgv {
	namespace gui {

		struct ellipsoid_render_ray_marching_style_gui_creator : public gui_creator
		{
			/// attempt to create a gui and return whether this was successful
			bool create(provider* p, const std::string& label,
				void* value_ptr, const std::string& value_type,
				const std::string& gui_type, const std::string& options, bool*)
			{
				if (value_type != cgv::type::info::type_name<cgv::render::ellipsoid_render_ray_marching_style>::get_name())
					return false;
				cgv::render::ellipsoid_render_ray_marching_style* rs_ptr = reinterpret_cast<cgv::render::ellipsoid_render_ray_marching_style*>(value_ptr);
				cgv::base::base* b = dynamic_cast<cgv::base::base*>(p);

				p->add_member_control(b, "Default Size", rs_ptr->size, "value_slider", "min=0.01;max=100;log=true;ticks=true");
				p->add_member_control(b, "Size Scale", rs_ptr->size_scale, "value_slider", "min=0.01;max=100;log=true;ticks=true");
				
				if (p->begin_tree_node("Cone Ray Marching Parameters", rs_ptr->rm, false)) {
					p->align("\a");
					p->add_member_control(b, "Epsilon", rs_ptr->rm.epsilon, "value_slider", "min=0;max=0.1;step=0.001;log=true");
					p->add_member_control(b, "Max Iterations", rs_ptr->rm.max_iterations, "value_slider", "min=0;max=50;step=1;unsigned=true");
					p->add_member_control(b, "FinDiffGrad Delta", rs_ptr->rm.fdg_delta, "value_slider", "min=0.01;max=1.0;step=0.01;ticks=true");
					p->add_member_control(b, "Show Bounding", rs_ptr->rm.show_bounding, "check");
					p->align("\b");
					p->end_tree_node(rs_ptr->rm);
				}
				return true;
			}
		};

#include <cgv_gl/gl/lib_begin.h>

		cgv::gui::gui_creator_registration<ellipsoid_render_ray_marching_style_gui_creator> ellipsoid_render_ray_marching_rs_gc_reg("ellipsoid_render_ray_marching_style_gui_creator");

	}
}