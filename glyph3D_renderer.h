#pragma once

#include "group_renderer.h"
//TODO THESIS
#include <cgv_reflect_types/media/illum/textured_surface_material.h>
#include <cgv/media/illum/textured_surface_material.h>

#include "gl/lib_begin.h"

namespace cgv { // @<
	namespace render { // @<

		/// color and opacity can be mapped independently to surface material
		enum GlyphColorMapping {
			GCM_NONE = 0,
			GCM_COLOR_FRONT = 1,
			GCM_COLOR_BACK = 2,
			GCM_COLOR = GCM_COLOR_FRONT | GCM_COLOR_BACK,
			GCM_OPACITY_FRONT = 4,
			GCM_OPACITY_BACK = 8,
			GCM_OPACITY = GCM_OPACITY_FRONT | GCM_OPACITY_BACK,
			GCM_COLOR_AND_OPACITY = GCM_COLOR | GCM_OPACITY
		};

		/** style of a point */
		struct /*CGV_API*/ glyph3D_render_style : public group_render_style
		{			
			/// default value for color when map color to material is used
			//TODO THESIS
			cgv::media::illum::surface_material::color_type surface_color;
			/// default value for the surface opacity when map color to material is used
			float surface_opacity;
			/// culling mode for point splats, set to CM_OFF in constructor
			CullingMode culling_mode;
			/// illumination mode defaults to \c IM_ONE_SIDED
			IlluminationMode illumination_mode;
			/// material side[s] where color is to be mapped to the diffuse material component, defaults to MS_FRONT_AND_BACK
			GlyphColorMapping map_color_to_material;
			/// material of surface
			//TODO THESIS
			cgv::media::illum::textured_surface_material material;
			///
			glyph3D_render_style();
		};

		/// base classes for renderers that support glyph3D rendering
		class /*CGV_API*/ glyph3D_renderer : public group_renderer
		{
		protected:
			//TODO THESIS
			bool has_normals;
			//TODO THESIS
			bool has_texcoords;
			bool cull_per_primitive;
		public:
			glyph3D_renderer();
			/// call this before setting attribute arrays to manage attribute array in given manager
			void enable_attribute_array_manager(const context& ctx, attribute_array_manager& aam);
			/// call this after last render/draw call to ensure that no other users of renderer change attribute arrays of given manager
			void disable_attribute_array_manager(const context& ctx, attribute_array_manager& aam);
			/// 
			bool enable(context& ctx);
			///
			bool disable(context& ctx);
			/// specify a single normal for all lines
			template <typename T>
			//TODO THESIS
			void set_normal(const context& ctx, const cgv::math::fvec<T, 3>& normal) { has_normals = true;  ref_prog().set_attribute(ctx, get_prog_attribute_location(ctx, "normal"), normal); }
			/// templated method to set the normal attribute array from a vector of normals of type T, which should have 3 components
			template <typename T>
			//TODO THESIS
			void set_normal_array(const context& ctx, const std::vector<T>& normals) { has_normals = true;  set_attribute_array(ctx, "normal", normals); }
			/// templated method to set the normal attribute from an array of normals of type T, which should have 3 components
			template <typename T>
			//TODO THESIS
			void set_normal_array(const context& ctx, const T* normals, size_t nr_elements, unsigned stride_in_bytes = 0) { has_normals = true;  set_attribute_array(ctx, "normal", normals, nr_elements, stride_in_bytes); }
			/// method to set the normal attribute from a vertex buffer object, the element type must be given as explicit template parameter
			//TODO THESIS
			void set_normal_array(const context& ctx, type_descriptor element_type, const vertex_buffer& vbo, size_t offset_in_bytes, size_t nr_elements, unsigned stride_in_bytes = 0);
			/// template method to set the normal attribute from a vertex buffer object, the element type must be given as explicit template parameter
			template <typename T>
			//TODO THESIS
			void set_normal_array(const context& ctx, const vertex_buffer& vbo, size_t offset_in_bytes, size_t nr_elements, unsigned stride_in_bytes = 0) { set_normal_array(ctx, type_descriptor(element_descriptor_traits<T>::get_type_descriptor(T()), true), vbo, offset_in_bytes, nr_elements, stride_in_bytes); }
			/// remove the normal attribute
			//TODO THESIS
			void remove_normal_array(const context& ctx);
			/// templated method to set the texcoord attribute without array
			template <typename T>
			//TODO THESIS
			void set_texcoord(const context& ctx, const T& texcoord) { has_texcoords = true;  ref_prog().set_attribute(ctx, get_prog_attribute_location(ctx, "texcoord"), texcoord); }
			/// templated method to set the texcoord attribute array from a vector of texcoords of type T
			template <typename T>
			//TODO THESIS
			void set_texcoord_array(const context& ctx, const std::vector<T>& texcoords) { has_texcoords = true;  set_attribute_array(ctx, "texcoord", texcoords); }
			/// templated method to set the texcoord attribute from an array of texcoords of type T
			template <typename T>
			//TODO THESIS
			void set_texcoord_array(const context& ctx, const T* texcoords, size_t nr_elements, unsigned stride_in_bytes = 0) { has_texcoords = true;  set_attribute_array(ctx, "texcoord", texcoords, nr_elements, stride_in_bytes); }
			/// template method to set the texcoord attribute from a vertex buffer object, the element type must be given as explicit template parameter
			//TODO THESIS
			void set_texcoord_array(const context& ctx, type_descriptor element_type, const vertex_buffer& vbo, size_t offset_in_bytes, size_t nr_elements, unsigned stride_in_bytes = 0);
			/// template method to set the texcoord attribute from a vertex buffer object, the element type must be given as explicit template parameter
			template <typename T>
			//TODO THESIS
			void set_texcoord_array(const context& ctx, const vertex_buffer& vbo, size_t offset_in_bytes, size_t nr_elements, unsigned stride_in_bytes = 0) { set_texcoord_array(ctx, type_descriptor(element_descriptor_traits<T>::get_type_descriptor(T()), true), vbo, offset_in_bytes, nr_elements, stride_in_bytes); }
			/// remove the texcoord attribute
			//TODO THESIS
			void remove_texcoord_array(const context& ctx);
		};
		struct /*CGV_API*/ glyph3D_render_style_reflect : public glyph3D_render_style
		{
			bool self_reflect(cgv::reflect::reflection_handler& rh);
		};
		extern /*CGV_API*/ cgv::reflect::extern_reflection_traits<glyph3D_render_style, glyph3D_render_style_reflect> get_reflection_traits(const glyph3D_render_style&);
	}
}

#include <cgv/config/lib_end.h>