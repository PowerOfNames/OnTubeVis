#pragma once

// C++ STL
#include <vector>
#include <set>

// CGV framework core
#include <cgv/base/node.h>
#include <cgv/gui/event_handler.h>
#include <cgv/gui/provider.h>
#include <cgv/render/drawable.h>

// CGV OpenGL lib
//#include <cgv_gl/rounded_cone_renderer.h>
//#include <cgv_gl/spline_tube_renderer.h>

// CGV framework graphics utility
#include <cgv_glutil/frame_buffer_container.h>
#include <cgv_glutil/radix_sort_4way.h>
#include <cgv_glutil/shader_library.h>

// local includes
#include "traj_loader.h"
#include "textured_spline_tube_renderer.h"

using namespace cgv::render;

////
// Plugin definition

/// baseline visualization plugin for arbitrary trajectory data as tubes using the framework tube renderers and
/// trajectory loading facilities
class tubes :
	public cgv::base::node,             // derive from node to integrate into global tree structure and to store a name
	public cgv::base::argument_handler, // derive from argument handler to be able to process custom arguments
	public cgv::gui::provider,          // derive from provider to obtain a GUI tab
	public cgv::gui::event_handler,     // derive from event handler to be able to directly react to user interaction
	public cgv::render::drawable        // derive from drawable for being able to render
{
public:
	/// real number type
	//typedef float real;

protected:
	/// store a pointer to the view for fast access
	view* view_ptr = nullptr;
	
	/// store the current OpenGL viewport configuration
	GLint viewport[4];

	/// path of the dataset to load - can be either a directory or a single file
	std::string datapath;

	/// misc configurable fields
	struct {
		/// proxy for controlling fltk_gl_view::instant_redraw
		bool instant_redraw_proxy = false;

		/// proxy for controlling stereo_view_interactor::fix_view_up_dir
		bool fix_view_up_dir_proxy = false;
	} misc_cfg;

	/// drag-n-drop state fields
	struct {
		/// current mouse position
		ivec2 pos;

		/// current drag-n-drop string
		std::string text;

		/// list of filenames extracted from @ref #text
		std::vector<std::string> filenames;
	} dnd;

	/// dataset state fields
	struct {
		/// set of filepaths for loading
		std::set<std::string> files;
	} dataset;

	cgv::glutil::frame_buffer_container fbc;
	cgv::glutil::shader_library shaders;

	/// rendering state fields
	struct {
		/// render style for the textured spline tubes
		textured_spline_tube_render_style style;
		
		/// accessor for the render data generated by the trajectory manager
		const traj_manager<float>::render_data *data;

		/// segment-wise arclength approximations (cubic bezier curves returning global
		/// trajectory arclength at the segment)
		std::vector<vec4> arclen_data;

		/// GPU-side storage buffer mirroring the \ref #arclen_data .
		vertex_buffer arclen_sbo;

		/// shared attribute array manager used by both renderers
		attribute_array_manager aam;

		/// the gpu sorter used to reorder the indices according to their corresponding segment visibility order
		cgv::glutil::gpu_sorter* sorter = nullptr;

		/// whether to sort the segemnts, which is used to boost performance together with conservative depth testing
		bool sort = true;

		float percentage = 1.0f;
	} render;

	/// trajectory manager
	traj_manager<float> traj_mgr;

	box3 bbox;
	
	/// test texture
	texture tex;
	texture density_tex;

	void set_view(void);
	void update_attribute_bindings(void);
	void calculate_bounding_box(void);

	float sd_quadratic_bezier(const vec3& A, const vec3& B, const vec3& C, const vec3& pos);
	std::vector<std::pair<int, float>> traverse_line(vec3& a, vec3& b, vec3& vbox_min, float vsize, ivec3& res);
	void create_density_volume(context& ctx, unsigned resolution);

	void set_ao_uniforms(context& ctx, const box3& volume_bbox, const uvec3 volume_resolution);

	/// draw methods
	void draw_dnd(context& ctx);
	void draw_trajectories(context& ctx);

public:
	tubes();
	std::string get_type_name() const { return "tubes"; }
	void handle_args(std::vector<std::string> &args);

	void clear(context& ctx);

	bool self_reflect(cgv::reflect::reflection_handler& rh);
	void stream_help(std::ostream& os);
	void stream_stats(std::ostream& os) {}

	bool handle(cgv::gui::event& e);
	void on_set(void* member_ptr);

	bool init(context& ctx);
	void init_frame(context& ctx);
	void draw(context& ctx);

	void create_gui();

	//bool on_exit_request();
};
