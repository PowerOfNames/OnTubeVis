#include <limits>
#include <cgv_gl/gl/gl.h>
#include <cgv_gl/gl/gl_tools.h>
#include "glyph3D_spline_tube_renderer.h"

namespace cgv {
	namespace render {
		glyph3D_spline_tube_renderer& ref_glyph3D_spline_tube_renderer(context& ctx, int ref_count_change, textured_spline_tube_render_style* textured_rs_ptr)
		{
			static int ref_count = 0;
			static glyph3D_spline_tube_renderer r;
			if(textured_rs_ptr != nullptr)
				r.set_textured_spline_tube_render_style_ptr(textured_rs_ptr);
			r.manage_singleton(ctx, "glyph3D_spline_tube_renderer", ref_count, ref_count_change);
			return r;
		};

		render_style* glyph3D_spline_tube_renderer::create_render_style() const
		{
			return new glyph3D_spline_tube_render_style();
		}

		glyph3D_spline_tube_render_style::glyph3D_spline_tube_render_style()
		{
			//THESIS2:
			glyph_method = GM_TUBE_SURFACE_PIPELINE;

			tube_backside_render.use_distance_check = true;
			tube_backside_render.clip_caps = true;
			tube_backside_render.distance_tolerance_factor = 0.00080f;

			cull_frontface = true;
			front_transparency = 0.2f;
			fresnel_reflection = 5.4f;
		}

		glyph3D_spline_tube_renderer::glyph3D_spline_tube_renderer()
		{
			has_node_ids = false;
			has_radii = false;
			has_tangents = false;
		}

		/// call this before setting attribute arrays to manage attribute array in given manager
		void glyph3D_spline_tube_renderer::enable_attribute_array_manager(const context& ctx, attribute_array_manager& aam)
		{
			surface_renderer::enable_attribute_array_manager(ctx, aam);
			if (has_attribute(ctx, "radius"))
				has_radii = true;
			if (has_attribute(ctx, "tangent"))
				has_tangents = true;
		}
		/// call this after last render/draw call to ensure that no other users of renderer change attribute arrays of given manager
		void glyph3D_spline_tube_renderer::disable_attribute_array_manager(const context& ctx, attribute_array_manager& aam)
		{
			surface_renderer::disable_attribute_array_manager(ctx, aam);
			has_radii = false;
			has_tangents = false;
		}
		bool glyph3D_spline_tube_renderer::validate_attributes(const context& ctx) const
		{
			// validate set attributes
			//bool res = surface_renderer::validate_attributes(ctx);
			//return res;

			if (!has_node_ids) {
				ctx.error("renderer::enable() node id attribute not set");
				return false;
			}
			return true;
		}
		void glyph3D_spline_tube_renderer::set_additional_defines(shader_define_map& defines) {
			additional_defines = defines;
		}
		void glyph3D_spline_tube_renderer::update_defines(shader_define_map& defines)
		{
			const textured_spline_tube_render_style* rs = textured_rs;
			const glyph3D_spline_tube_render_style& rs3D = get_style<glyph3D_spline_tube_render_style>();

			defines.clear();

			//THESIS: (Danke David)
			//shader_code::set_define(defines, "GLYPH_TYPE_IS_3D", rs.glyph_dimension, rs.GD_2D);
			//THESIS2:
			shader_code::set_define(defines, "BACKSIDE_DEPTH", rs3D.glyph_method, rs3D.GM_TUBE_SURFACE_PIPELINE);

			shader_code::set_define(defines, "USE_CONSERVATIVE_DEPTH", rs->use_conservative_depth, false);
			if (rs->is_tube()) {
				shader_code::set_define(defines, "USE_CUBIC_TANGENTS", rs->use_cubic_tangents, true);
				shader_code::set_define(defines, "USE_VIEW_SPACE_POSITION", rs->use_view_space_position, true);
				shader_code::set_define(defines, "PRIMITIVE_INTERSECTOR", rs->line_primitive, rs->LP_TUBE_RUSSIG);
				//THESIS2:
				shader_code::set_define(defines, "BACKSIDE_DEPTH", rs3D.glyph_method, rs3D.GM_TUBE_SURFACE_PIPELINE);
				static const bool no = false;
				shader_code::set_define(defines, "USE_RIBBONS", no, false);

				//THESIS:
				// glyph type
				/*const bool is3D = rs3D.is_3D();
				shader_code::set_define(defines, "GLYPH_TYPE_IS_3D", is3D, false);*/
				//THESIS2:
				const bool isBB = rs3D.is_billboard_pipeline();
				shader_code::set_define(defines, "BACKSIDE_DEPTH", isBB, false);				
			}
			else if (rs->line_primitive == rs->LP_RIBBON_GEOMETRY) {
				static const bool yes = true;
				shader_code::set_define(defines, "USE_RIBBONS", yes, false);
			}
			shader_code::set_define(defines, "ATTRIB_MODE", rs->attrib_mode, rs->AM_ALL);
			shader_code::set_define(defines, "MODE", rs->fragment_mode, rs->FM_RAY_CAST);
			if (rs->line_primitive != rs->LP_RIBBON_GEOMETRY)
				shader_code::set_define(defines, "BOUNDING_GEOMETRY_TYPE", rs->bounding_geometry, rs->BG_ALIGNED_BOX_BILLBOARD);
			if (rs->line_primitive == rs->LP_RIBBON_RAYCASTED) {
				shader_code::set_define(defines, "EXACT_RIBBON_BBOXES", rs->rcribbon.exact_ribbon_bboxes, false);
				shader_code::set_define(defines, "BBOX_COORD_SYSTEM", rs->rcribbon.bbox_coord_system, rs->rcribbon.BBO_RCC);
				shader_code::set_define(defines, "RAY_CENTRIC_ISECTS", rs->rcribbon.ray_centric_isects, false);
				shader_code::set_define(defines, "MAX_INTERSECTION_STACK_SIZE", rs->rcribbon.max_intersection_stack_size, (unsigned)8);
				shader_code::set_define(defines, "DBG_VISUALIZE_STATS", rs->rcribbon.debug.visualize_stats, rs->rcribbon.debug.VS_OFF);
				shader_code::set_define(defines, "DBG_VISUALIZE_LEAF_BBOXES", rs->rcribbon.debug.visualize_leaf_bboxes, false);
			}

			

			for (const auto& define : additional_defines)
				defines.insert(define);
		}
		bool glyph3D_spline_tube_renderer::build_shader_program(context& ctx, shader_program& prog, const shader_define_map& defines)
		{
			const textured_spline_tube_render_style* rs = textured_rs;
			const glyph3D_spline_tube_render_style& rs3D = get_style<glyph3D_spline_tube_render_style>();
			last_active_glyph_method = rs3D.glyph_method;

			//THESIS:
			//TODO THESIS:
			if (rs->is_tube())
			{
				return prog.build_program(ctx, "spline_tube_glyph3D.glpr", true, defines);
			}
			//Not implemented in Thesis
			else if (rs->line_primitive == rs->LP_RIBBON_RAYCASTED)
				return prog.build_program(ctx, "view_aligned_ribbon.glpr", true, defines);
			else
				return prog.build_program(ctx, "textured_spline_ribbon.glpr", true, defines);
			
		}
		bool glyph3D_spline_tube_renderer::enable(context& ctx)
		{
			const textured_spline_tube_render_style* rs = textured_rs;
			const glyph3D_spline_tube_render_style& rs3D = get_style<glyph3D_spline_tube_render_style>();

			//THESIS2: needed?
			if (last_active_glyph_method != rs3D.glyph_method) {
				clear(ctx);
				init(ctx);
			}

			if (!surface_renderer::enable(ctx))
				return false;

			if (!ref_prog().is_linked())
				return false;

			ref_prog().set_uniform(ctx, "radius_scale", rs->radius_scale);
			ref_prog().set_uniform(ctx, "cyclopic_eye", cyclopic_eye);
			ref_prog().set_uniform(ctx, "view_dir", view_dir);
			ref_prog().set_uniform(ctx, "viewport", viewport);
			ref_prog().set_uniform(ctx, "cap_clip_distance", rs->cap_clip_distance);
			ref_prog().set_uniform(ctx, "max_t", rs->max_t);

			if (rs3D.is_billboard_pipeline())
			{
				ref_prog().set_uniform(ctx, "tbr.use_distance_check", rs3D.tube_backside_render.use_distance_check);
				ref_prog().set_uniform(ctx, "tbr.clip_caps", rs3D.tube_backside_render.clip_caps);
				ref_prog().set_uniform(ctx, "tbr.distance_tolerance_factor", rs3D.tube_backside_render.distance_tolerance_factor);
				ref_prog().set_uniform(ctx, "cull_frontface", rs3D.cull_frontface);
			}

			if (rs->line_primitive == rs->LP_RIBBON_RAYCASTED) {
				ref_prog().set_uniform(ctx, "linearity_thr", rs->rcribbon.linearity_thr);
				ref_prog().set_uniform(ctx, "screwiness_thr", std::min(rs->rcribbon.screwiness_thr, .9921875f));
				ref_prog().set_uniform(ctx, "subdiv_abort_thr", rs->rcribbon.subdiv_abort_thr);
			}

			return true;
		}
		///
		bool glyph3D_spline_tube_renderer::disable(context& ctx)
		{
			if (!attributes_persist()) {
				has_radii = false;
				has_tangents = false;
			}

			return surface_renderer::disable(ctx);
		}

		void glyph3D_spline_tube_renderer::draw(context& ctx, size_t start, size_t count, bool use_strips, bool use_adjacency, uint32_t strip_restart_index)
		{
			glDisable(GL_CULL_FACE);
			draw_impl(ctx, PT_POINTS, start, count, use_strips, use_adjacency, strip_restart_index);
			glEnable(GL_CULL_FACE);
		}

		bool glyph3D_spline_tube_render_style_reflect::self_reflect(cgv::reflect::reflection_handler& rh)
		{
			return
				rh.reflect_base(*static_cast<surface_render_style*>(this)) &&
				rh.reflect_member("glyph_method", glyph_method) &&
				rh.reflect_member("tube_backside_render_cull_frontface", cull_frontface) &&
				rh.reflect_member("tube_backside_render_use_distance_check", tube_backside_render.use_distance_check) &&
				rh.reflect_member("tube_backside_render_clip_caps", tube_backside_render.clip_caps) &&
				rh.reflect_member("tube_backside_render_distance_tolerance_factor", tube_backside_render.distance_tolerance_factor) &&
				rh.reflect_member("front_transparency", front_transparency) &&
				rh.reflect_member("fresnel_reflection", fresnel_reflection);
		}

		cgv::reflect::extern_reflection_traits<glyph3D_spline_tube_render_style, glyph3D_spline_tube_render_style_reflect> get_reflection_traits(const glyph3D_spline_tube_render_style&)
		{
			return cgv::reflect::extern_reflection_traits<glyph3D_spline_tube_render_style, glyph3D_spline_tube_render_style_reflect>();
		}
	}
}

namespace cgv {
	namespace reflect {	

		//THESIS2
		enum_reflection_traits<cgv::render::glyph3D_spline_tube_render_style::GlyphMethod> get_reflection_traits(const cgv::render::glyph3D_spline_tube_render_style::GlyphMethod&) {
			return enum_reflection_traits<cgv::render::glyph3D_spline_tube_render_style::GlyphMethod>("GM_TUBE_SURFACE_PIPELINE,GM_BILLBOARD_PIPELINE");
		}
	}
}

#include <cgv/gui/provider.h>

namespace cgv {
	namespace gui {

		struct glyph3D_spline_tube_render_style_gui_creator : public gui_creator {
			/// attempt to create a gui and return whether this was successful
			bool create(provider* p, const std::string& label,
				void* value_ptr, const std::string& value_type,
				const std::string& gui_type, const std::string& options, bool*) {
				if (value_type != cgv::type::info::type_name<cgv::render::glyph3D_spline_tube_render_style>::get_name())
					return false;
				cgv::render::glyph3D_spline_tube_render_style* rs_ptr = reinterpret_cast<cgv::render::glyph3D_spline_tube_render_style*>(value_ptr);
				cgv::base::base* b = dynamic_cast<cgv::base::base*>(p);

				if (rs_ptr->is_billboard_pipeline())
				{
					if (p->begin_tree_node("Tube Backside Parameters", rs_ptr->tube_backside_render, false)) {
						p->align("\a");
						p->add_member_control(b, "Cull Frontface", rs_ptr->cull_frontface, "check");
						p->add_member_control(b, "Use Distance Checks", rs_ptr->tube_backside_render.use_distance_check, "check");
						p->add_member_control(b, "Clip Caps", rs_ptr->tube_backside_render.clip_caps, "check");
						p->add_member_control(b, "Distance Tolerance", rs_ptr->tube_backside_render.distance_tolerance_factor, "value_slider", "min=0;max=0.1;step=0.00001;log=true");
						p->align("\b");
						p->end_tree_node(rs_ptr->tube_backside_render);
					}
					p->add_member_control(b, "Front Transparency", rs_ptr->front_transparency, "value_slider", "min=0.0;max=1.0;step=0.01;ticks=true");
					p->add_member_control(b, "Fresnel RefracIdx", rs_ptr->fresnel_reflection, "value_slider", "min=1.0;max=10.0;step=0.1;ticks=true");
				}
				return true;
			}
		};

#include <cgv_gl/gl/lib_begin.h>

		cgv::gui::gui_creator_registration<glyph3D_spline_tube_render_style_gui_creator> glyph3D_spline_tube_rs_gc_reg("glyph3D_spline_tube_render_style_gui_creator");
	}
}
