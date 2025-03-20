#include "cone_renderer_ray_marching.h"


#include <cgv_gl/gl/gl.h>
#include <cgv_gl/gl/gl_tools.h>

namespace cgv {
	namespace render {
		cone_renderer_ray_marching& ref_cone_renderer_ray_marching(context& ctx, int ref_count_change)
		{
			static int ref_count = 0;
			static cone_renderer_ray_marching r;
			r.manage_singleton(ctx, "cone_renderer_ray_marching", ref_count, ref_count_change);
			return r;
		}

		render_style* cone_renderer_ray_marching::create_render_style() const
		{
			return new cone_render_ray_marching_style();
		}

		cone_render_ray_marching_style::cone_render_ray_marching_style()
		{
			radius = 1.0f;
			radius_scale = 1.0f;

			show_caps = true;
			composite_arrow = false;
		}

		cone_renderer_ray_marching::cone_renderer_ray_marching()
		{
			has_radii = false;
			shader_defines = shader_define_map();
		}

		/// call this before setting attribute arrays to manage attribute array in given manager
		void cone_renderer_ray_marching::enable_attribute_array_manager(const context& ctx, attribute_array_manager& aam)
		{
			surface_renderer::enable_attribute_array_manager(ctx, aam);
			if (has_attribute(ctx, "radius"))
				has_radii = true;
		}
		/// call this after last render/draw call to ensure that no other users of renderer change attribute arrays of given manager
		void cone_renderer_ray_marching::disable_attribute_array_manager(const context& ctx, attribute_array_manager& aam)
		{
			surface_renderer::disable_attribute_array_manager(ctx, aam);
			has_radii = false;
		}
		void cone_renderer_ray_marching::remove_radius_array(const context& ctx) {
			has_radii = false;
			remove_attribute_array(ctx, "radius");
		}
		bool cone_renderer_ray_marching::validate_attributes(const context& ctx) const
		{
			const cone_render_ray_marching_style& crs = get_style<cone_render_ray_marching_style>();
			bool res = surface_renderer::validate_attributes(ctx);
			return res;
		}
		void cone_renderer_ray_marching::update_defines(shader_define_map& defines)
		{
			const cone_render_ray_marching_style& crs = get_style<cone_render_ray_marching_style>();
			shader_code::set_define(defines, "CAPS", crs.show_caps, true);	
			shader_code::set_define(defines, "COMP_ARROW", crs.composite_arrow, false);		
		}
		bool cone_renderer_ray_marching::build_shader_program(context& ctx, shader_program& prog, const shader_define_map& defines)
		{
			return prog.build_program(ctx, "cone_ray_marched.glpr", true, defines);
		}		
		/// 
		bool cone_renderer_ray_marching::enable(context& ctx)
		{
			if (!surface_renderer::enable(ctx))
				return false;

			if (!ref_prog().is_linked())
				return false;

			const cone_render_ray_marching_style& crs = get_style<cone_render_ray_marching_style>();
			if (!has_radii)
				ref_prog().set_attribute(ctx, "radius", crs.radius);

			ref_prog().set_uniform(ctx, "radius_scale", crs.radius_scale);
			
			return true;
		}
		///
		bool cone_renderer_ray_marching::disable(context& ctx)
		{
			
			if (!attributes_persist()) {
				has_radii = false;
			}
			return surface_renderer::disable(ctx);
		}

		bool cone_render_ray_marching_style_reflect::self_reflect(cgv::reflect::reflection_handler& rh)
		{
			return
				rh.reflect_base(*static_cast<surface_render_style*>(this)) &&
				rh.reflect_member("radius", radius);
		}

		void cone_renderer_ray_marching::draw(context& ctx, size_t start, size_t count, bool use_strips, bool use_adjacency, uint32_t strip_restart_index)
		{
			draw_impl(ctx, PT_LINES, start, count, use_strips, use_adjacency, strip_restart_index);
		}

		void cone_renderer_ray_marching::clear(const context& ctx)
		{
			renderer::clear(ctx);
		}

		cgv::reflect::extern_reflection_traits<cone_render_ray_marching_style, cone_render_ray_marching_style_reflect> get_reflection_traits(const cone_render_ray_marching_style&)
		{
			return cgv::reflect::extern_reflection_traits<cone_render_ray_marching_style, cone_render_ray_marching_style_reflect>();
		}
	}
}

#include <cgv/gui/provider.h>

namespace cgv {
	namespace gui {

		struct cone_render_ray_marching_style_gui_creator : public gui_creator {
			/// attempt to create a gui and return whether this was successful
			bool create(provider* p, const std::string& label,
				void* value_ptr, const std::string& value_type,
				const std::string& gui_type, const std::string& options, bool*) {
				if (value_type != cgv::type::info::type_name<cgv::render::cone_render_ray_marching_style>::get_name())
					return false;
				cgv::render::cone_render_ray_marching_style* crs_ptr = reinterpret_cast<cgv::render::cone_render_ray_marching_style*>(value_ptr);
				cgv::base::base* b = dynamic_cast<cgv::base::base*>(p);

				p->add_member_control(b, "Default Radius", crs_ptr->radius, "value_slider", "min=0.001;step=0.0001;max=10.0;log=true;ticks=true");
				p->add_member_control(b, "Radius Scale", crs_ptr->radius_scale, "value_slider", "min=0.01;step=0.0001;max=100.0;log=true;ticks=true");

				p->add_member_control(b, "Show Caps", crs_ptr->show_caps, "check");
				p->add_member_control(b, "Composite Arrow", crs_ptr->composite_arrow, "check");
				
				p->add_gui("surface_render_style", *static_cast<cgv::render::surface_render_style*>(crs_ptr));
				return true;
			}
		};

#include <cgv_gl/gl/lib_begin.h>

		cgv::gui::gui_creator_registration<cone_render_ray_marching_style_gui_creator> cone_rm_rs_gc_reg("cone_render_ray_marching_style_gui_creator");
	}
}
