/*
 * PyMOL COLLADA export
 *
 * (c) 2014 Jared Sampson
 * (c) 2014 Schrodinger, Inc.
 */

#include "os_predef.h"
#include "os_std.h"

#include "Feedback.h"
#include "Ray.h"

#ifdef _HAVE_LIBXML

#include <libxml/xmlwriter.h>
#include <libxml/encoding.h>

#include <time.h>

#include "Setting.h"
#include "Ortho.h"
#include "Color.h"
#include "Util.h"
#include "Version.h"
#include "MemoryDebug.h"
#include "Sphere.h"
#include "Scene.h"
#include "CGO.h"


#define XML_VERSION "1.0"
#define XML_ENCODING "UTF-8"

#define PRECISION 0.001f
#define TRANS_PRECISION 0.01f

/*
 * Size of individual text nodes must be less than the maximum size for
 * the LIBXML Parser node: 10MB. Text nodes longer than this will be split.
 */
#define XML_NODE_SIZE_LIMIT 1000000

/* For debugging output. */
#define COLLADA_DEBUG 0

enum {
  COLLADA_GEOM_MODE_DEFAULT = 0,
  COLLADA_GEOM_MODE_BLENDER = 1,

  /* Total number of `collada_geometry_mode` settings implemented. */
  COLLADA_GEOM_MODE_COUNT
};

/* Returns a standard XML-formatted timestamp */
static char *GetTimestamp(){
  time_t now;
  struct tm *local;
  char buffer[20];

  time(&now);
  local = localtime(&now);
  strftime(buffer, 20, "%Y-%m-%dT%H:%M:%S", local);
  return strdup(buffer);
}

/*
 * Returns index of the first float in an array matching the given value to
 * the given precision.
 */
static int GetFloatPositionInArray(float val, float *arr, int len,
    float precision)
{
  int i;
  for (i = 0; i < len; i++) {
    if (precision > fabsf(val - arr[i])){
      return i;
    }
  }
  return -1;  // not found
}


/* Returns 1 if any of the 4 VLA strings has gotten too big, else 0. */
static int ReachedXmlNodeSizeLimit(int pos_str_cc, int norm_str_cc,
    int col_str_cc, int p_str_cc)
{
  int r = 0;
  if(pos_str_cc >= XML_NODE_SIZE_LIMIT)
    r = 1;
  if(norm_str_cc >= XML_NODE_SIZE_LIMIT)
    r = 1;
  if(col_str_cc >= XML_NODE_SIZE_LIMIT)
    r = 1;
  if(p_str_cc >= XML_NODE_SIZE_LIMIT)
    r = 1;
  return r;
}


/* Writes global <asset> element with file meta information. */
static void ColladaWriteAssetElement(xmlTextWriterPtr w, int geom_mode)
{
  xmlTextWriterStartElement(w, BAD_CAST "asset");
  xmlTextWriterStartElement(w, BAD_CAST "contributor");
  xmlTextWriterWriteElement(w, BAD_CAST "author",
      BAD_CAST "PyMOL User");
#ifdef _PyMOL_VERSION
  char *vers = (char *)"PyMOL " _PyMOL_VERSION;
#else
  char *vers = (char *)"PyMOL";
#endif
  xmlTextWriterWriteElement(w, BAD_CAST "authoring_tool",
      BAD_CAST vers);
  xmlTextWriterEndElement(w);  // contributor

  char *ts = GetTimestamp();

  xmlTextWriterWriteElement(w, BAD_CAST "created", BAD_CAST ts);
  xmlTextWriterWriteElement(w, BAD_CAST "modified", BAD_CAST ts);
  free(ts);

  if (geom_mode != COLLADA_GEOM_MODE_BLENDER) {
    xmlTextWriterStartElement(w, BAD_CAST "unit");
    xmlTextWriterWriteAttribute(w, BAD_CAST "meter", BAD_CAST "1e-10");
    xmlTextWriterWriteAttribute(w, BAD_CAST "name", BAD_CAST "Angstrom");
    xmlTextWriterEndElement(w);  // unit
  }

  xmlTextWriterWriteElement(w, BAD_CAST "up_axis", BAD_CAST "Y_UP");
  xmlTextWriterEndElement(w);  // asset
}


/* Writes the <library_cameras> element, including PyMOL's viewport camera. */
static void ColladaWriteLibraryCameras(xmlTextWriterPtr w, PyMOLGlobals *G, int width,
    int height, float fov, float front, float back)
{
  SceneViewType view;
  SceneGetView(G, view);

  float aspect_ratio = (float)width / (float)height;
  int ortho = SettingGetGlobal_i(G, cSetting_ortho);
  int ray_ortho = SettingGetGlobal_i(G, cSetting_ray_orthoscopic);

  if (ray_ortho == -1) {
    ray_ortho = ortho;
  }

  xmlTextWriterStartElement(w, BAD_CAST "library_cameras");

  xmlTextWriterStartElement(w, BAD_CAST "camera");
  xmlTextWriterWriteAttribute(w, BAD_CAST "id", BAD_CAST "camera");
  xmlTextWriterStartElement(w, BAD_CAST "optics");
  xmlTextWriterStartElement(w, BAD_CAST "technique_common");

  if (ray_ortho) {
    float ymag = SettingGetGlobal_i(G, cSetting_field_of_view) / 2;
    float ymag_scale = -view[18] / 50;  // scale fov based on camera distance
    float ortho_adjust = 0.88f;  // prevent apparent size change
    ymag = ymag * ymag_scale * ortho_adjust;

    xmlTextWriterStartElement(w, BAD_CAST "orthographic");
    xmlTextWriterWriteFormatElement(w, BAD_CAST "ymag", "%6.4f", ymag);
  } else {
    xmlTextWriterStartElement(w, BAD_CAST "perspective");
    xmlTextWriterWriteFormatElement(w, BAD_CAST "yfov", "%6.4f", fov);
  }

  xmlTextWriterWriteFormatElement(w, BAD_CAST "aspect_ratio",
      "%6.4f", aspect_ratio);
  xmlTextWriterWriteFormatElement(w, BAD_CAST "znear", "%6.4f", front);
  xmlTextWriterWriteFormatElement(w, BAD_CAST "zfar", "%6.4f", back);
  xmlTextWriterEndElement(w);  // perspective or orthographic

  xmlTextWriterEndElement(w);  // technique_common
  xmlTextWriterEndElement(w);  // optics
  xmlTextWriterEndElement(w);  // camera

  xmlTextWriterEndElement(w);  // library_cameras
}


/* Writes the <library_lights> element, including ambient and directional
 * PyMOL lights. */
static void ColladaWriteLibraryLights(xmlTextWriterPtr w, CRay *I, PyMOLGlobals *G)
{
  xmlTextWriterStartElement(w, BAD_CAST "library_lights");

  /* Ambient light */
  float ambient = SettingGetGlobal_f(I->G, cSetting_ambient);
  if(ambient > 0.5){
    ambient = 0.5;
  }

  /* RGB = "0.5 0.5 0.5" seems to work well; `ambient` setting is handled in
   * default effect */
  xmlTextWriterStartElement(w, BAD_CAST "light");
  xmlTextWriterWriteAttribute(w, BAD_CAST "id", BAD_CAST "ambient-light");
  xmlTextWriterStartElement(w, BAD_CAST "technique_common");
  xmlTextWriterStartElement(w, BAD_CAST "ambient");
  xmlTextWriterWriteElement(w, BAD_CAST "color", BAD_CAST "0.5 0.5 0.5");
  xmlTextWriterEndElement(w);  // directional
  xmlTextWriterEndElement(w);  // technique_common
  xmlTextWriterEndElement(w);  // light

  /* All other lights are directional, and their intensities are equal. Use
   * instances of this light in the <node> section below. */
  xmlTextWriterStartElement(w, BAD_CAST "light");
  xmlTextWriterWriteAttribute(w, BAD_CAST "id", BAD_CAST "pymol-light");
  xmlTextWriterStartElement(w, BAD_CAST "technique_common");
  xmlTextWriterStartElement(w, BAD_CAST "directional");

  /* Intensity decreases with increasing number of lights. */
  float intensity = SceneGetSpecularValue(G, 0.6, 10);
  char *value = (char *) malloc(sizeof(char) * 30);
  sprintf(value, "%5.3f %5.3f %5.3f", intensity, intensity, intensity);
  xmlTextWriterWriteElement(w, BAD_CAST "color", BAD_CAST value);
  free(value);

  xmlTextWriterEndElement(w);  // directional
  xmlTextWriterEndElement(w);  // technique_common
  xmlTextWriterEndElement(w);  // light

  xmlTextWriterEndElement(w);  // library_lights
}


/* Writes a <color>R G B</color> element. */
static void ColladaWriteCommonColorElement(xmlTextWriterPtr w, char *name,
    char *sid, char *value)
{
  xmlTextWriterStartElement(w, BAD_CAST name);
  xmlTextWriterStartElement(w, BAD_CAST "color");
  if(sid) {
    xmlTextWriterWriteAttribute(w, BAD_CAST "sid", BAD_CAST sid);
  }
  else {
    xmlTextWriterWriteAttribute(w, BAD_CAST "sid", BAD_CAST name);
  }
  xmlTextWriterWriteString(w, BAD_CAST value);
  xmlTextWriterEndElement(w);  // color
  xmlTextWriterEndElement(w);  // element
}


/* Writes a <float>value</float> element. */
static void ColladaWriteCommonFloatElement(xmlTextWriterPtr w, char *name,
    char *sid, char *value)
{
  xmlTextWriterStartElement(w, BAD_CAST name);
  xmlTextWriterStartElement(w, BAD_CAST "float");
  if(sid) {
    xmlTextWriterWriteAttribute(w, BAD_CAST "sid", BAD_CAST sid);
  }
  else {
    xmlTextWriterWriteAttribute(w, BAD_CAST "sid", BAD_CAST name);
  }
  xmlTextWriterWriteString(w, BAD_CAST value);
  xmlTextWriterEndElement(w);  // float
  xmlTextWriterEndElement(w);  // element
}


/* Writes a set of <phong> shader parameter elements. */
static void ColladaWritePhongEffect(xmlTextWriterPtr w, char *id,
    float amb, float spec, float shin,
    float trans, float iref)
{
  char *value = (char *) malloc(100 * sizeof(char));

  xmlTextWriterStartElement(w, BAD_CAST "effect");
  xmlTextWriterWriteAttribute(w, BAD_CAST "id",
      BAD_CAST id);
  xmlTextWriterStartElement(w, BAD_CAST "profile_COMMON");
  xmlTextWriterStartElement(w, BAD_CAST "technique");
  xmlTextWriterWriteAttribute(w, BAD_CAST "sid",
      BAD_CAST "common");

  xmlTextWriterStartElement(w, BAD_CAST "phong");

  if (amb > PRECISION) {
    sprintf(value, "0.5 0.5 0.5 %5.3f", amb);
    ColladaWriteCommonColorElement(w, (char *)"ambient", NULL, value);
  }

  if (spec > PRECISION) {
    sprintf(value, "0.5 0.5 0.5 %5.3f", spec);
    ColladaWriteCommonColorElement(w, (char *)"specular", NULL, value);
  }

  if(shin > PRECISION) {
    sprintf(value, "%5.3f", shin);
    ColladaWriteCommonFloatElement(w, (char *)"shininess", NULL, value);
  }

  if(trans > PRECISION) {
    sprintf(value, "%5.3f", trans);
    ColladaWriteCommonFloatElement(w, (char *)"transparency", NULL, value);
  }

  if(iref > PRECISION) {
    sprintf(value, "%5.3f", iref);
    ColladaWriteCommonFloatElement(w, (char *)"index_of_refraction", NULL, value);
  }

  xmlTextWriterEndElement(w);  // phong
  xmlTextWriterEndElement(w);  // technique
  xmlTextWriterEndElement(w);  // profile_COMMON
  xmlTextWriterEndElement(w);  // effect

  free(value);
}


/* Writes a <material> element as an instance of a specific <effect> by URL. */
static void ColladaWriteInstanceEffectMaterial(xmlTextWriterPtr w,
    char *id, char *url)
{
  xmlTextWriterStartElement(w, BAD_CAST "material");
  xmlTextWriterWriteAttribute(w, BAD_CAST "id", BAD_CAST id);
  xmlTextWriterStartElement(w, BAD_CAST "instance_effect");
  xmlTextWriterWriteAttribute(w, BAD_CAST "url", BAD_CAST url);
  xmlTextWriterEndElement(w);  // instance_effect
  xmlTextWriterEndElement(w);  // material
}


/* Opens a <geometry id="..."><mesh> element. */
static void ColladaBeginGeometryMesh(xmlTextWriterPtr w, int geom){
  xmlTextWriterStartElement(w, BAD_CAST "geometry");
  xmlTextWriterWriteFormatAttribute(w, BAD_CAST "id", "geom%i", geom);
  xmlTextWriterStartElement(w, BAD_CAST "mesh");
}


/* Closes a </mesh></geometry> element. */
static void ColladaEndGeometryMesh(xmlTextWriterPtr w){
  xmlTextWriterEndElement(w);  // mesh
  xmlTextWriterEndElement(w);  // geometry
}


/* Writes the <library_effects> element, including transparency values. */
static void ColladaWriteLibraryEffects(xmlTextWriterPtr w, PyMOLGlobals *G,
    int trans_len, float *trans)
{
  xmlTextWriterStartElement(w, BAD_CAST "library_effects");

  float amb = SettingGetGlobal_f(G, cSetting_ambient);
  if(amb > 0.5){
    amb = 0.5;
  }

  float spec = SettingGetGlobal_f(G, cSetting_specular_intensity);

  float shin_factor = 5.0f;
  float shin = SettingGetGlobal_f(G, cSetting_shininess) / shin_factor;

  /* Default effect */
  ColladaWritePhongEffect(w, (char *)"default-effect", amb, spec, shin, 1, 1);

  /* Background effect */
  ColladaWritePhongEffect(w, (char *)"bg-effect", 0.5f, 0, 0, 0, 0);

  /* Transparency effects */
  int i;
  char *name = (char *) malloc(100 * sizeof(char));
  for (i = 0; i < trans_len; i++) {
    sprintf(name, "transparency-%1.2f-effect", trans[i]);
    ColladaWritePhongEffect(w, name, amb, spec, shin, 1 - trans[i], 1);
#if COLLADA_DEBUG
    printf("Wrote Phong effect: %s\n", name);
#endif
  }
  xmlTextWriterEndElement(w);  // library_effects
  free(name);
}


/* Writes the <library_materials> element, including transparent materials. */
static void ColladaWriteLibraryMaterials(xmlTextWriterPtr w, int trans_len, float *trans)
{
  xmlTextWriterStartElement(w, BAD_CAST "library_materials");

  /* Default material */
  ColladaWriteInstanceEffectMaterial(w, (char *)"default-material",
      (char *)"#default-effect");
  ColladaWriteInstanceEffectMaterial(w, (char *)"bg-material",
      (char *)"#bg-effect");

  /* Transparent materials */
  int i;
  char *name = (char *) malloc(100 * sizeof(char));
  char *url = (char *) malloc(100 * sizeof(char));
  for (i = 0; i < trans_len; i++) {

    sprintf(name, "transparency-%1.2f-material", trans[i]);
    sprintf(url, "#transparency-%1.2f-effect", trans[i]);

    ColladaWriteInstanceEffectMaterial(w, name, url);
#if COLLADA_DEBUG
    printf("Wrote material: %s for effect: %s\n", name, url);
#endif

  }
  xmlTextWriterEndElement(w);  // library_materials
  free(name);
  free(url);
}


/*
 * Writes a 3-dimensional <source> element, e.g. XYZ or RGB, to define
 * position, normal, or color source data.
 *
 * Note: "count" should be the number of accessible elements (e.g. vertices,
 * normals; NOT individual float values).  The length of the float array will
 * be 3x this value.
 */
static void ColladaWrite3DSource(xmlTextWriterPtr w, int geom, char *name, int count,
    char *data_str, char *dim)
{
  xmlTextWriterStartElement(w, BAD_CAST "source");
  xmlTextWriterWriteFormatAttribute(w, BAD_CAST "id",
      "geom%i-mesh-%s", geom, name);

  xmlTextWriterStartElement(w, BAD_CAST "float_array");
  xmlTextWriterWriteFormatAttribute(w, BAD_CAST "id",
      "geom%i-mesh-%s-array", geom, name);
  xmlTextWriterWriteFormatAttribute(w, BAD_CAST "count", "%i", count * 3);
  xmlTextWriterWriteString(w, BAD_CAST data_str);
  xmlTextWriterEndElement(w);  // float_array

  xmlTextWriterStartElement(w, BAD_CAST "technique_common");
  xmlTextWriterStartElement(w, BAD_CAST "accessor");
  xmlTextWriterWriteFormatAttribute(w, BAD_CAST "source",
      "#geom%i-mesh-%s-array", geom, name);
  xmlTextWriterWriteFormatAttribute(w, BAD_CAST "count", "%i", count);
  xmlTextWriterWriteAttribute(w, BAD_CAST "stride", BAD_CAST "3");

  char n[2] = "\0";
  sprintf(n, "%c", dim[0]);
  xmlTextWriterStartElement(w, BAD_CAST "param");
  xmlTextWriterWriteAttribute(w, BAD_CAST "name", BAD_CAST n);
  xmlTextWriterWriteAttribute(w, BAD_CAST "type", BAD_CAST "float");
  xmlTextWriterEndElement(w);

  sprintf(n, "%c", dim[1]);
  xmlTextWriterStartElement(w, BAD_CAST "param");
  xmlTextWriterWriteAttribute(w, BAD_CAST "name", BAD_CAST n);
  xmlTextWriterWriteAttribute(w, BAD_CAST "type", BAD_CAST "float");
  xmlTextWriterEndElement(w);

  sprintf(n, "%c", dim[2]);
  xmlTextWriterStartElement(w, BAD_CAST "param");
  xmlTextWriterWriteAttribute(w, BAD_CAST "name", BAD_CAST n);
  xmlTextWriterWriteAttribute(w, BAD_CAST "type", BAD_CAST "float");
  xmlTextWriterEndElement(w);

  xmlTextWriterEndElement(w);  // accessor
  xmlTextWriterEndElement(w);  // technique_common
  xmlTextWriterEndElement(w);  // source
}


/* Write a <vertices> element using the current geometry's positions source. */
static void ColladaWriteVertices(xmlTextWriterPtr w, int geom)
{
  xmlTextWriterStartElement(w, BAD_CAST "vertices");
  xmlTextWriterWriteFormatAttribute(w, BAD_CAST "id",
      "geom%i-mesh-vertices", geom);
  xmlTextWriterStartElement(w, BAD_CAST "input");
  xmlTextWriterWriteAttribute(w, BAD_CAST "semantic",
      BAD_CAST "POSITION");
  xmlTextWriterWriteFormatAttribute(w, BAD_CAST "source",
      "#geom%i-mesh-positions", geom);
  xmlTextWriterEndElement(w);  // input
  xmlTextWriterEndElement(w);  // vertices
}


/* Writes vertex, normal, and color <input> elements for the given geometry. */
static void ColladaWriteVNCInputs(xmlTextWriterPtr w, int geom)
{
  xmlTextWriterStartElement(w, BAD_CAST "input");
  xmlTextWriterWriteAttribute(w, BAD_CAST "offset", BAD_CAST "0");
  xmlTextWriterWriteAttribute(w, BAD_CAST "semantic",
      BAD_CAST "VERTEX");
  xmlTextWriterWriteFormatAttribute(w, BAD_CAST "source",
      "#geom%i-mesh-vertices", geom);
  xmlTextWriterEndElement(w);  // input

  xmlTextWriterStartElement(w, BAD_CAST "input");
  xmlTextWriterWriteAttribute(w, BAD_CAST "offset", BAD_CAST "1");
  xmlTextWriterWriteAttribute(w, BAD_CAST "semantic",
      BAD_CAST "NORMAL");
  xmlTextWriterWriteFormatAttribute(w, BAD_CAST "source",
      "#geom%i-mesh-normals", geom);
  xmlTextWriterEndElement(w);  // input

  xmlTextWriterStartElement(w, BAD_CAST "input");
  xmlTextWriterWriteAttribute(w, BAD_CAST "offset", BAD_CAST "2");
  xmlTextWriterWriteAttribute(w, BAD_CAST "semantic",
      BAD_CAST "COLOR");
  xmlTextWriterWriteFormatAttribute(w, BAD_CAST "source",
      "#geom%i-mesh-colors", geom);
  xmlTextWriterEndElement(w);  // input
}


/* Writes a <vcount> element with the given string as its value. */
static void ColladaWriteVCountElement(xmlTextWriterPtr w, char *vcount_str)
{
  xmlTextWriterStartElement(w, BAD_CAST "vcount");
  xmlTextWriterWriteFormatString(w, "%s", vcount_str);
  xmlTextWriterEndElement(w);
}


/* Writes a <vcount> element for the given number of triangles. */
static void ColladaWriteTrianglesVCountElement(xmlTextWriterPtr w, int tri)
{
  int i;
  char *vc_str = VLACalloc(char, 1000);
  ov_size cc = 0;
  char *next = (char *) malloc(10 * sizeof(char));

  for(i = 0; i < tri; i++){
    sprintf(next, "3 ");  // all triangles
    UtilConcatVLA(&vc_str, &cc, next);
  }

  ColladaWriteVCountElement(w, vc_str);

  VLAFree(vc_str);
  free(next);
}


/* Writes a <p> element with the given string as its value. */
static void ColladaWritePrimitiveElement(xmlTextWriterPtr w, char *p_str)
{
  xmlTextWriterStartElement(w, BAD_CAST "p");
  xmlTextWriterWriteFormatString(w, "%s", p_str);
  xmlTextWriterEndElement(w);  // p
}


/* Writes a <triangles> element with the current primitive (<p>) string. */
static void ColladaWriteTrianglesElement(xmlTextWriterPtr w, int geom, int tri,
    char *p_str, int mode = COLLADA_GEOM_MODE_DEFAULT)
{
  xmlTextWriterStartElement(w,
      BAD_CAST(mode == COLLADA_GEOM_MODE_BLENDER ? "polylist" : "triangles"));

  xmlTextWriterWriteFormatAttribute(w, BAD_CAST "count", "%i", tri);
  xmlTextWriterWriteFormatAttribute(w, BAD_CAST "material",
      "geom%i-material", geom);

  ColladaWriteVNCInputs(w, geom);

  if (mode == COLLADA_GEOM_MODE_BLENDER) {
    ColladaWriteTrianglesVCountElement(w, tri);
  }

  ColladaWritePrimitiveElement(w, p_str);

  xmlTextWriterEndElement(w);  // triangles
}


/* Writes a <polylist> element consisting entirely of triangles. */
static void ColladaWriteTrianglesPolylistElement(xmlTextWriterPtr w, int geom, int tri,
    char *p_str)
{
  ColladaWriteTrianglesElement(w, geom, tri, p_str, COLLADA_GEOM_MODE_BLENDER);
}


/*
 * Opens a <polylist> element, including <input>s.  Must be followed by at
 * least one <p> element and then closed.
 */
static void ColladaBeginPolylistElement(xmlTextWriterPtr w, int geom, int count)
{
  xmlTextWriterStartElement(w, BAD_CAST "polylist");
  xmlTextWriterWriteFormatAttribute(w, BAD_CAST "count", "%i", count);
  xmlTextWriterWriteFormatAttribute(w, BAD_CAST "material",
      "geom%i-material", geom);

  ColladaWriteVNCInputs(w, geom);
}


/* Closes a <polylist> element. */
static void ColladaEndPolylistElement(xmlTextWriterPtr w)
{
  xmlTextWriterEndElement(w);
}


/*
 * Opens a <tristrips> element, including inputs.  Must be followed by at
 * least one <p> element.  Multiple <p> elements can be used to have
 * multiple triangle strips using the same vertices.
 */
static void ColladaBeginTristripsElement(xmlTextWriterPtr w, int geom, int num_strips)
{
  xmlTextWriterStartElement(w, BAD_CAST "tristrips");
  xmlTextWriterWriteFormatAttribute(w, BAD_CAST "count",
      "%i", num_strips);
  xmlTextWriterWriteFormatAttribute(w, BAD_CAST "material",
      "geom%i-material", geom);

  ColladaWriteVNCInputs(w, geom);
}


/* Closes a <tristrips> element. */
static void ColladaEndTristripsElement(xmlTextWriterPtr w)
{
  xmlTextWriterEndElement(w);
}


#if 0
/*
 * Opens a <trifans> element, including inputs.  Must be followed by at
 * least one <p> element.
 */
static void ColladaBeginTrifansElement(xmlTextWriterPtr w, int geom, int num_fans)
{
  xmlTextWriterStartElement(w, BAD_CAST "trifans");
  xmlTextWriterWriteFormatAttribute(w, BAD_CAST "count",
      "%i", num_fans);
  xmlTextWriterWriteFormatAttribute(w, BAD_CAST "material",
      "geom%i-material", geom);

  ColladaWriteVNCInputs(w, geom);
}


/* Closes a <trifans> element. */
static void ColladaEndTrifansElement(xmlTextWriterPtr w)
{
  xmlTextWriterEndElement(w);
}
#endif


/* Writes a complete <geometry> element for a triangle mesh. */
static void ColladaWriteMeshGeometry(xmlTextWriterPtr w, int geom,
    int pos, char* positions_str,
    int norm, char* normals_str,
    int col, char* colors_str,
    int tri, char* p_str, int mode)
{
  ColladaBeginGeometryMesh(w, geom);
  ColladaWrite3DSource(w, geom, (char *)"positions", pos, positions_str,
      (char *)"XYZ");
  ColladaWrite3DSource(w, geom, (char *)"normals", norm, normals_str,
      (char *)"XYZ");
  ColladaWrite3DSource(w, geom, (char *)"colors", col, colors_str,
      (char *)"RGB");
  ColladaWriteVertices(w, geom);

  ColladaWriteTrianglesElement(w, geom, tri, p_str, mode);

  ColladaEndGeometryMesh(w);
}

#endif // _HAVE_LIBXML

/* Generates COLLADA output and appends it to `vla_ptr`. */
void RayRenderCOLLADA(CRay * I, int width, int height,
    char **vla_ptr, float front, float back,
    float fov)
{
#ifndef _HAVE_LIBXML
  PRINTFB(I->G, FB_Ray, FB_Errors)
    " ColladaRender-Error: No libxml2 support, can't output COLLADA file.\n"
    ENDFB(I->G);
#elif !(defined(LIBXML_WRITER_ENABLED) && defined(LIBXML_OUTPUT_ENABLED))
  PRINTFB(I->G, FB_Ray, FB_Errors)
    " ColladaRender-Error: libxml2 library not properly configured " \
    "with writer and output support, can't output COLLADA file.\n"
    ENDFB(I->G);
#else

  char *vla = *vla_ptr;
  ov_size cc = 0;               /* character count */
  PyMOLGlobals *G = I->G;


  /*
   * Setting: geometry_export_mode
   * -----------------------------
   * 0 = Output everything. (default)
   * 1 = Output geometry and materials only; exclude lighting and camera
   *     information from scene.
   */
  int identity = (SettingGetGlobal_i(I->G, cSetting_geometry_export_mode) == 1);

  /*
   * Setting: collada_background_box
   * -------------------------------
   * 0 = Do not include a background box. (default)
   * 1 = Include the background box for more accurate reproduction of the scene.
   * NB: This setting is overridden by geometry_export_mode == 1, in which case
   * no background box will be generated.
   */
  int bgbox = SettingGetGlobal_i(I->G, cSetting_collada_background_box);

  /*
   * Setting: collada_geometry_mode
   * ------------------------------
   * 0 = Valid COLLADA 1.4.1. (default)
   * 1 = Blender-compatible (only <polylist> geometries are used).
   */
  int geom_mode = SettingGetGlobal_i(I->G, cSetting_collada_geometry_mode);
  if (geom_mode >= COLLADA_GEOM_MODE_COUNT || geom_mode < 0) {
    geom_mode = COLLADA_GEOM_MODE_DEFAULT;
  }

  /*
   * Setting: collada_export_lighting
   * --------------------------------
   * 0 = No lighting information included in output.  Best for
   *     interactive scenes. (default)
   * 1 = Lighting information (<library_lights> and light <node>s) included
   *     in output.
   */
  int lighting = (SettingGetGlobal_i(I->G, cSetting_collada_export_lighting) == 1);
  int lc = SettingGetGlobal_i(I->G, cSetting_light_count);


  /* Ray trace */
  RayExpandPrimitives(I);
  RayTransformFirst(I, 0, identity);
  RayComputeBox(I);

  /* initialize for XML writing */
  int rc;  // return codes for error handling
  xmlTextWriterPtr w;
  xmlDocPtr doc;
  xmlChar *tmp;
  int buffersize;

  /* Create a new XML DOM tree, to which the COLLADA document will be
   * written */
  doc = xmlNewDoc(BAD_CAST XML_VERSION);
  if (doc == NULL) {
    printf("ColladaRender: Error creating the xml document tree (xmlNewDoc).\n");
    return;
  }

  /* Create a new XmlWriter */
  w = xmlNewTextWriterTree(doc, NULL, 0);
  if (w == NULL) {
    printf("ColladaRender: Error creating the xml writer (xmlNewTextWriterTree).\n");
    return;
  }

  /* Start the XML document */
  rc = xmlTextWriterStartDocument(w, NULL, XML_ENCODING, NULL);
  if (rc < 0) {
    printf("ColladaRender: Error at xmlTextWriterStartDocument\n");
    return;
  }

  /* Begin COLLADA */
  xmlTextWriterStartElement(w, BAD_CAST "COLLADA");
  xmlTextWriterWriteAttribute(w, BAD_CAST "xmlns",
      BAD_CAST "http://www.collada.org/2005/11/COLLADASchema");
  xmlTextWriterWriteAttribute(w, BAD_CAST "version", BAD_CAST "1.4.1");

  /* Asset */
  ColladaWriteAssetElement(w, geom_mode);

  /* Cameras */
  if (!identity) {
    ColladaWriteLibraryCameras(w, G, width, height, fov, front, back);
  }

  /* Lights */
  if (!identity && lighting) {
    ColladaWriteLibraryLights(w, I, G);
  }

  /* Geometries */
  int geom = 0;

  /*
   * Track transparency levels used on a per-geom basis to be added to
   * additional library_effects and library_materials elements.
   */
  float *trans = (float *) malloc(sizeof(float));  // store transparency values
  int trans_len = 0;

  /* Associate geometries with transparency values. */
  int *geom_trans = (int *) malloc(sizeof(int));
  int geom_trans_len = 1;

  {
    xmlTextWriterStartElement(w, BAD_CAST "library_geometries");

    int a, tri, pos, norm, col;
    CPrimitive *prim;
    int mesh_obj = false;
    int largest_dim = 10;
    float cur_trans = 0.0f;

    /* Initialize data string VLAs and character counts. */
    char *positions_str = VLACalloc(char, 1000);
    char   *normals_str = VLACalloc(char, 1000);
    char    *colors_str = VLACalloc(char, 1000);
    char         *p_str = VLACalloc(char, 1000);
    ov_size  pos_str_cc = 0;
    ov_size norm_str_cc = 0;
    ov_size  col_str_cc = 0;
    ov_size    p_str_cc = 0;

    /*
     * Loop through primitives, plus one extra time to finish writing a
     * triangle mesh if necessary.
     */
    for (a = 0; a <= I->NPrimitive; a++) {
      prim = I->Primitive + a;

      /* Handle transitions between/after triangle meshes. */
      if (mesh_obj) {
        if (a == I->NPrimitive ||
            prim->type != cPrimTriangle ||
            TRANS_PRECISION <= fabsf(prim->trans - cur_trans))
        {

          /*
           * Not a triangle primitive, but the previous mesh is still
           * active; OR, transparency level is different; OR the final
           * primitive was part of a triangle mesh that needs to be
           * closed.
           *
           * Write the previous triangle mesh from its data strings
           * and counters, reset mesh tracking variables.
           */

          ColladaWriteMeshGeometry(w, geom, pos, positions_str,
              norm, normals_str, col, colors_str, tri, p_str, geom_mode);

          geom += 1;
          mesh_obj = false;
        }
      }
      if(!mesh_obj) {
        /* First triangle primitive or any other primititve. */

        /* Allocate data strings */
        pos_str_cc = 0;
        norm_str_cc = 0;
        col_str_cc = 0;
        p_str_cc = 0;


        if(a < I->NPrimitive){

          /* Reset counters */
          tri = 0;
          pos = 0;
          norm = 0;
          col = 0;
          cur_trans = prim->trans;

          /* Increase geom_trans with geom, realloc exponentially. */
          while (geom >= geom_trans_len) {
            geom_trans_len *= 2;
            geom_trans = (int *) realloc(geom_trans,
                (geom_trans_len) * sizeof(int));
#if COLLADA_DEBUG
            printf("  increased geom_trans_len to %i\n", geom_trans_len);
#endif
          }

          /* Record transparency for each geometry. */
          int p = GetFloatPositionInArray(cur_trans, trans,
              trans_len, TRANS_PRECISION);

          if (0 > p) {
            /* Not found, add new transparency level. */
            trans = (float *) realloc(trans, (trans_len + 1) * sizeof(float));
            trans[trans_len] = cur_trans;
#if COLLADA_DEBUG
            printf("  Added new trans[%i] = %1.2f\n", trans_len,
                trans[trans_len]);
#endif
            geom_trans[geom] = trans_len;
            trans_len++;
          } else {
            /* Reference existing transparency level. */
            geom_trans[geom] = p;
          }

#if COLLADA_DEBUG
          printf("geom_trans[%i] = %i\n", geom, geom_trans[geom]);
#endif

          if(prim->type == cPrimTriangle){

            /* Track the open triangle mesh. */
            mesh_obj = true;
          }
        }
        else {
          /*
           * No active mesh object, and the last primitive has already been
           * processed.
           */

          /* Generate bounding box geometry for background. */
          /* Only if collada_background_box == 1 and geometry_export_mode == 0 */
          if (!identity && bgbox) {
            /* Ensure camera is inside bounding box */
            int i;
            for (i = 0; i < 3; i++) {
              if (fabsf(I->Pos[i]) > largest_dim) {
                largest_dim = fabsf(I->Pos[i]);
              }
            }

            /* Leave plenty of room */
            largest_dim *= 100;

            /* Allocate data strings */
            pos_str_cc = 0;
            norm_str_cc = 0;
            col_str_cc = 0;
            p_str_cc = 0;

            int cube_coords[24] = {
              1,  1,  1,
              1, -1,  1,
              -1, -1,  1,
              -1,  1,  1,
              -1,  1, -1,
              -1, -1, -1,
              1, -1, -1,
              1,  1, -1 };

            int cube_face_verts[24] = {
              0, 1, 2, 3,
              3, 2, 5, 4,
              4, 5, 6, 7,
              7, 6, 1, 0,
              0, 3, 4, 7,
              1, 6, 5, 2 };

            char *next = (char *) malloc(200 * sizeof(int));

            /* positions */
            for(i = 0; i < 24; i++){
              sprintf(next, "%i ", cube_coords[i] * largest_dim);
              UtilConcatVLA(&positions_str, &pos_str_cc, next);
            }

            /* normals */
            for(i = 0; i < 24; i++){
              sprintf(next, "%i ", -cube_coords[i]);
              UtilConcatVLA(&normals_str, &norm_str_cc, next);
            }

            /* color */
            const float *bg_color;
            bg_color = ColorGet(G, SettingGet_color(G, NULL, NULL, cSetting_bg_rgb));

            sprintf(next, "%6.4f %6.4f %6.4f", bg_color[0], bg_color[1], bg_color[2]);
            UtilConcatVLA(&colors_str, &col_str_cc, next);

            /* p */
            for(i = 0; i < 6; i++){
              sprintf(next, "%i %i 0 %i %i 0 %i %i 0 %i %i 0  ",
                  cube_face_verts[4 * i], cube_face_verts[4 * i],
                  cube_face_verts[4 * i + 1], cube_face_verts[4 * i + 1],
                  cube_face_verts[4 * i + 2], cube_face_verts[4 * i + 2],
                  cube_face_verts[4 * i + 3], cube_face_verts[4 * i + 3]);
              UtilConcatVLA(&p_str, &p_str_cc, next);
            }

            xmlTextWriterStartElement(w, BAD_CAST "geometry");
            xmlTextWriterWriteAttribute(w, BAD_CAST "id", BAD_CAST "geom-bg");
            xmlTextWriterStartElement(w, BAD_CAST "mesh");

            ColladaWrite3DSource(w, geom, (char *)"positions", 8,
                positions_str, (char *)"XYZ");
            ColladaWrite3DSource(w, geom, (char *)"normals", 8,
                normals_str, (char *)"XYZ");
            ColladaWrite3DSource(w, geom, (char *)"colors", 1,
                colors_str, (char *)"RGB");
            ColladaWriteVertices(w, geom);

            xmlTextWriterStartElement(w, BAD_CAST "polylist");
            xmlTextWriterWriteAttribute(w, BAD_CAST "count", BAD_CAST "6");
            xmlTextWriterWriteAttribute(w, BAD_CAST "material", BAD_CAST "geom-bg-material");

            ColladaWriteVNCInputs(w, geom);
            xmlTextWriterWriteElement(w, BAD_CAST "vcount", BAD_CAST "4 4 4 4 4 4");
            ColladaWritePrimitiveElement(w, p_str);
            xmlTextWriterEndElement(w);  // polylist

            xmlTextWriterEndElement(w);  // mesh
            xmlTextWriterEndElement(w);  // geometry

            free(next);

          }


          /* Finished with geometries. */
          break;
        }
      }

      /* Update largest dimension - assumes distance from origin for most
       * objects will be much greater than e.g. sphere/cylinder/cone radius
       * and camera won't be "too far" away. */
      {
        int i;
        for(i = 0; i < 3; i++) {
          switch (prim->type) {
            /* 3 vertices defined */
            case cPrimTriangle:
              if (largest_dim < prim->v3[i]) {
                largest_dim = prim->v3[i];
              }
              /* 2 vertices defined */
            case cPrimCone:
            case cPrimCylinder:
              if (largest_dim < prim->v2[i]) {
                largest_dim = prim->v2[i];
              }
              /* only 1 vertex defined */
            case cPrimSphere:
            default:
              if (largest_dim < prim->v1[i]) {
                largest_dim = prim->v1[i];
              }
              break;
          }
        }
      }

      /* Process the primitive according to its type. */
      switch(prim->type){
        case cPrimSphere:         // 1
          {
#if COLLADA_DEBUG > 1
            printf("Primitive %i: cPrimSphere\n", a);
#endif

            /* sphere_quality range: 0 to (NUMBER_OF_SPHERE_LEVELS - 1). */
            int sq = (SettingGetGlobal_i(I->G, cSetting_sphere_quality));

            SphereRec *sp;
            int i, j;
            int *q, *s;
            char *next = (char *) malloc(200 * sizeof(char));  // enough for 9 %6.4 floats

            sp = I->G->Sphere->Sphere[sq];
            q = sp->Sequence;
            s = sp->StripLen;


            /* Generate <source> elements */

            /* For each tristrip, add positions and normals to their respective
             * data strings. */
            for(i = 0; i < sp->NStrip; i++) {
#if COLLADA_DEBUG > 1
              printf("\ni = %i, StripLen = %i\n", i, *s);
#endif

              for(j = 0; j < (*s); j++) {
                /* positions */
                sprintf(next, "%6.4f %6.4f %6.4f ",
                    prim->v1[0] + (prim->r1 * sp->dot[*q][0]),
                    prim->v1[1] + (prim->r1 * sp->dot[*q][1]),
                    prim->v1[2] + (prim->r1 * sp->dot[*q][2]));
#if COLLADA_DEBUG > 1
                printf("positions j = %02i: %s\n", j, next);
#endif
                UtilConcatVLA(&positions_str, &pos_str_cc, next);

                /* normals */
                sprintf(next, "%6.4f %6.4f %6.4f ",
                    sp->dot[*q][0], sp->dot[*q][1], sp->dot[*q][2]);
                UtilConcatVLA(&normals_str, &norm_str_cc, next);
                q++;  // next position in sequence

              }
              s++;  // next strip

            }

            /* Colors: only one color per sphere. */
            sprintf(next, "%6.4f %6.4f %6.4f",
                prim->c1[0], prim->c1[1], prim->c1[2]);
            UtilConcatVLA(&colors_str, &col_str_cc, next);



            ColladaBeginGeometryMesh(w, geom);
            ColladaWrite3DSource(w, geom, (char *)"positions", sp->NVertTot,
                positions_str, (char *)"XYZ");
            ColladaWrite3DSource(w, geom, (char *)"normals", sp->NVertTot,
                normals_str, (char *)"XYZ");
            ColladaWrite3DSource(w, geom, (char *)"colors", 1,
                colors_str, (char *)"RGB");
            ColladaWriteVertices(w, geom);

            /* Reset to the beginning of the list of strips. */
            q = sp->Sequence;
            s = sp->StripLen;

            int vert_ct;
            vert_ct = 0;

            if (geom_mode == COLLADA_GEOM_MODE_BLENDER) {

              /* polylists */
              int tri = 0;

              for (i = 0; i < sp->NStrip; i++) {

                /* Each triangle has the last 2 vertices of previous one as its
                 * first 2 vertices.  Color is the same for all vertices. */
                for(j = 2; j < (*s); j++) {

                  if (j % 2 == 0) {
                    /* even, normals correct as-is */
                    sprintf(next, "%i %i 0 %i %i 0 %i %i 0 ",
                        vert_ct + j - 2, vert_ct + j - 2,
                        vert_ct + j - 1, vert_ct + j - 1,
                        vert_ct + j, vert_ct + j);
                  } else {
                    /* odd, switch order to maintain correct normal */
                    sprintf(next, "%i %i 0 %i %i 0 %i %i 0 ",
                        vert_ct + j - 1, vert_ct + j - 1,
                        vert_ct + j - 2, vert_ct + j - 2,
                        vert_ct + j, vert_ct + j);
                  }

                  UtilConcatVLA(&p_str, &p_str_cc, next);
                  q++;
                  tri++;
                }

                ColladaWriteTrianglesPolylistElement(w, geom, tri, p_str);

                vert_ct += (*s);
                s++;

                /* Get a fresh p_str for the next strip. */
                VLAFree(p_str);
                p_str = VLACalloc(char, 1000);
                p_str_cc = 0;
                tri = 0;

              }


            } else {

              /* tristrips */
              ColladaBeginTristripsElement(w, geom, sp->NStrip);

              /* For each strip, write a <p> element. */
              for(i = 0; i < sp->NStrip; i++) {

                /* Write the 2 initial vertices. */
                sprintf(next, "%i %i 0 %i %i 0 ",
                    vert_ct, vert_ct, vert_ct + 1, vert_ct + 1);
                UtilConcatVLA(&p_str, &p_str_cc, next);

                /* After the first 2, each new vertex adds a new triangle.
                 * Color is the same for all vertices. */
                for(j = 2; j < (*s); j++) {
                  sprintf(next, "%i %i 0 ", vert_ct + j, vert_ct + j);
                  UtilConcatVLA(&p_str, &p_str_cc, next);
                  q++;
                }

                ColladaWritePrimitiveElement(w, p_str);

                /* Start the next strip at the next vertex. */
                vert_ct += *s;
                s++;

                /* Get a fresh p_str for the next strip. */
                VLAFree(p_str);
                p_str = VLACalloc(char, 1000);
                p_str_cc = 0;

              }

              ColladaEndTristripsElement(w);
            }

            free(next);

            ColladaEndGeometryMesh(w);
            geom += 1;
            break;

          }  // cPrimSphere

        case cPrimCylinder:       // 2
        case cPrimSausage:        // 4
        case cPrimCone:           // 7
          {
#if COLLADA_DEBUG > 1
            if(prim->type == cPrimCylinder) {
              printf("Primitive %i: cPrimCylinder\n", a);
            } else if(prim->type == cPrimSausage) {
              printf("Primitive %i: cPrimSausage\n", a);
            } else {
              printf("Primitive %i: cPrimCone\n", a);
            }
#endif

#define DAE_MAX_EDGE 50

            /* Local vars */
            float d[3], t[3], p0[3], p1[3], p2[3];
            float v_buf[9], *v, vv1[3], vv2[3], vvv1[3], vvv2[3];
            float x[DAE_MAX_EDGE + 1], y[DAE_MAX_EDGE + 1];
            float overlap, overlap2, nub, nub2, r2 = 0.F;
            int nEdge, c; //colorFlag;
            int i;
            char *next = (char *) malloc(200 * sizeof(char));

            char *cap1_p_str = VLACalloc(char, 1000);
            ov_size cap1_cc = 0;

            bool stick_round_nub = SettingGetGlobal_i(I->G, cSetting_stick_round_nub);
            const int j_arr[] = {2, 3, 2};
            int nCapTri = 0;
            cCylCap captype[2] = {prim->cap1, prim->cap2};
            CGO *cgocap[2] = {NULL, NULL};

            if(prim->type == cPrimSausage) {
              captype[0] = captype[1] = cCylCapRound;
            }

            v = v_buf;
            nEdge = SettingGetGlobal_i(I->G, cSetting_stick_quality);
            overlap = prim->r1 * SettingGetGlobal_f(I->G, cSetting_stick_overlap);
            nub = prim->r1 * SettingGetGlobal_f(I->G, cSetting_stick_nub);
            if(prim->type == cPrimCone) {
              r2 = prim->r2;
              overlap2 = prim->r2 * SettingGetGlobal_f(I->G, cSetting_stick_overlap);
              nub2 = prim->r2 * SettingGetGlobal_f(I->G, cSetting_stick_nub);
            }
            else {
              r2 = prim->r1;
              overlap2 = overlap;
              nub2 = nub;
            }

            if(nEdge > DAE_MAX_EDGE)
              nEdge = DAE_MAX_EDGE;
            subdivide(nEdge, x, y);

            //colorFlag = (prim->c1 != prim->c2) && prim->c2;

            /* primary axis vector p0 */
            p0[0] = (prim->v2[0] - prim->v1[0]);
            p0[1] = (prim->v2[1] - prim->v1[1]);
            p0[2] = (prim->v2[2] - prim->v1[2]);

            normalize3f(p0);

            /* cap1 */
            copy3f(prim->v1, vv1);
            copy3f(vv1, vvv1);

            /* cap2 */
            copy3f(prim->v2, vv2);
            copy3f(vv2, vvv2);

            d[0] = (vv2[0] - vv1[0]);
            d[1] = (vv2[1] - vv1[1]);
            d[2] = (vv2[2] - vv1[2]);

            get_divergent3f(d, t);

            /* orthogonal vectors p1, p2 */
            cross_product3f(d, t, p1);
            normalize3f(p1);

            cross_product3f(d, p1, p2);
            normalize3f(p2);

            /* now we have a coordinate system */

            /* cap1 */
            if(captype[0] == cCylCapRound) {
              if (stick_round_nub) {
                cgocap[0] = CGONew(I->G);
                CGORoundNub(cgocap[0], prim->v1, p0, p1, p2, -1, nEdge, prim->r1);
              } else {
                for(i = 0; i < 3; i++) {
                  vv1[i] -= p0[i] * overlap;
                  vvv1[i] = vv1[i] - p0[i] * nub;
                }
              }
            }

            /* cap2 */
            if(captype[1] == cCylCapRound) {
              if (stick_round_nub) {
                cgocap[1] = CGONew(I->G);
                CGORoundNub(cgocap[1], prim->v2, p0, p1, p2, 1, nEdge, prim->r1);
              } else {
                for(i = 0; i < 3; i++) {
                  vv2[i] += p0[i] * overlap2;
                  vvv2[i] = vv2[i] + p0[i] * nub2;
                }
              }
            }

            /* colors */
            sprintf(next, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f ",
                prim->c1[0], prim->c1[1], prim->c1[2],
                prim->c2[0], prim->c2[1], prim->c2[2]);
            UtilConcatVLA(&colors_str, &col_str_cc, next);

            /* Generate the data strings */
            char *vcount_str = VLACalloc(char, 100);
            ov_size vc_cc = 0;
            {

              /* Write values for each edge. The first edge is also the last
               * edge. */
              for(c = nEdge; c >= 0; c--) {

                /* vector only, not positioned yet */
                v[0] = p1[0] * x[c] + p2[0] * y[c];
                v[1] = p1[1] * x[c] + p2[1] * y[c];
                v[2] = p1[2] * x[c] + p2[2] * y[c];

                /* vertices */
                v[3] = vv1[0] + v[0] * prim->r1;
                v[4] = vv1[1] + v[1] * prim->r1;
                v[5] = vv1[2] + v[2] * prim->r1;
                v[6] = vv2[0] + v[0] * r2;
                v[7] = vv2[1] + v[1] * r2;
                v[8] = vv2[2] + v[2] * r2;


                /* positions */
                sprintf(next, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f ",
                    v[3], v[4], v[5], v[6], v[7], v[8]);
                UtilConcatVLA(&positions_str, &pos_str_cc, next);

                /* normals: Both vertices have the same normal. */
                sprintf(next, "%6.4f %6.4f %6.4f ",
                    v[0], v[1], v[2]);
                UtilConcatVLA(&normals_str, &norm_str_cc, next);

                if (c > 0) {
                  /* vcount */
                  sprintf(next, "4 ");
                  UtilConcatVLA(&vcount_str, &vc_cc, next);

                  /* p for polylist */
                  sprintf(next, "%i %i %i %i %i %i %i %i %i %i %i %i ",
                      pos, norm, col,
                      pos + 1, norm, col + 1,
                      pos + 3, norm + 1, col,
                      pos + 2, norm + 1, col + 1);
                  UtilConcatVLA(&p_str, &p_str_cc, next);
                }
                pos += 2;
                norm += 1;

              }

              /* Caps */
              {
                /* add another vertex-normal-color set for the center of each
                 * cap if r > 0 */

#define CONE_MIN_RADIUS 10e-6f

                if(prim->r1 > CONE_MIN_RADIUS && !cgocap[0] && captype[0] != cCylCapNone){
                  /* positions */
                  sprintf(next, "%6.4f %6.4f %6.4f ",
                      vvv1[0], vvv1[1], vvv1[2]);
                  UtilConcatVLA(&positions_str, &pos_str_cc, next);

                  /* normals */
                  sprintf(next, "%6.4f %6.4f %6.4f ",
                      -p0[0], -p0[1], -p0[2]);
                  UtilConcatVLA(&normals_str, &norm_str_cc, next);

                  /* p */
                  for(i = 0; i < nEdge; i++) {
                    sprintf(next, "%i %i %i %i %i %i %i %i %i ",
                        pos, norm, col,
                        2*i, i, col,
                        2*i+2, i+1, col);
                    UtilConcatVLA(&cap1_p_str, &cap1_cc, next);
                    ++nCapTri;
                  }

                  ++pos;
                  ++norm;
                }

                if(r2 > CONE_MIN_RADIUS && !cgocap[1] && captype[1] != cCylCapNone){
                  /* positions */
                  sprintf(next, "%6.4f %6.4f %6.4f ",
                      vvv2[0], vvv2[1], vvv2[2]);
                  UtilConcatVLA(&positions_str, &pos_str_cc, next);

                  /* normals */
                  sprintf(next, "%6.4f %6.4f %6.4f ",
                      p0[0], p0[1], p0[2]);
                  UtilConcatVLA(&normals_str, &norm_str_cc, next);

                  /* p */
                  for(i = 0; i < nEdge; i++) {
                    /* reverse order for other end */
                    sprintf(next, "%i %i %i %i %i %i %i %i %i ",
                        pos, norm, col + 1,
                        2*i+3, i+1, col + 1,
                        2*i+1, i, col + 1);
                    UtilConcatVLA(&cap1_p_str, &cap1_cc, next);
                    ++nCapTri;
                  }

                  ++pos;
                  ++norm;
                }
              }

#if COLLADA_DEBUG > 1
              printf("positions: %s\n", positions_str);
              printf("normals:   %s\n", normals_str);
              printf("colors:    %s\n", colors_str);
              printf("p:         %s\n", p_str);
#endif

              for (i = 0; i < 2; ++i) {
                if (!cgocap[i])
                  continue;

                int pos0 = -1;

                for (auto it = cgocap[i]->begin(); !it.is_stop(); ++it) {
                  auto pc = it.data();
                  int op = it.op_code();
                  switch (op){
                    case CGO_BEGIN:
                      pos0 = pos;
                      break;
                    case CGO_NORMAL:
                      /* normals */
                      sprintf(next, "%6.4f %6.4f %6.4f ", pc[0], pc[1], pc[2]);
                      UtilConcatVLA(&normals_str, &norm_str_cc, next);
                      ++norm;
                      break;
                    case CGO_VERTEX:
                      /* positions */
                      sprintf(next, "%6.4f %6.4f %6.4f ", pc[0], pc[1], pc[2]);
                      UtilConcatVLA(&positions_str, &pos_str_cc, next);
                      ++pos;

                      /* unroll the triangle strip */
                      if (pos > pos0 + 2) {
                        const int * j = j_arr + ((pos - pos0) % 2);
                        sprintf(next, "%i %i %i %i %i %i %i %i %i ",
                            pos - j[0], norm - j[0], col + i,
                            pos - j[1], norm - j[1], col + i,
                            pos - 1,    norm - 1,    col + i);
                        UtilConcatVLA(&cap1_p_str, &cap1_cc, next);
                        ++nCapTri;
                      }

                      break;
                  }
                }

                CGOFree(cgocap[i]);
              }
            }

            col += 2;

            ColladaBeginGeometryMesh(w, geom);
            ColladaWrite3DSource(w, geom, (char *)"positions", pos,
                positions_str, (char *)"XYZ");
            ColladaWrite3DSource(w, geom, (char *)"normals", norm,
                normals_str, (char *)"XYZ");
            ColladaWrite3DSource(w, geom, (char *)"colors", col,
                colors_str, (char *)"RGB");
            ColladaWriteVertices(w, geom);

            /* polylist for cylinder shaft */
            ColladaBeginPolylistElement(w, geom, nEdge);
            ColladaWriteVCountElement(w, vcount_str);
            ColladaWritePrimitiveElement(w, p_str);
            ColladaEndPolylistElement(w);

            if (nCapTri) {
              ColladaWriteTrianglesElement(w, geom, nCapTri, cap1_p_str, geom_mode);
            }

            ColladaEndGeometryMesh(w);
            geom += 1;


            VLAFree(vcount_str);
            VLAFree(cap1_p_str);

            free(next);

            break;
          }  // cPrimCylinder

        case cPrimTriangle:       // 3
          {
#if COLLADA_DEBUG > 1
            printf("Primitive %i: cPrimTriangle\n", a);
#endif

            char *next = (char *) malloc(200 * sizeof(char));  // enough for 9 color floats

            /*** Positions ***/
            sprintf(next, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f ",
                prim->v1[0], prim->v1[1], prim->v1[2],
                prim->v2[0], prim->v2[1], prim->v2[2],
                prim->v3[0], prim->v3[1], prim->v3[2]);
            UtilConcatVLA(&positions_str, &pos_str_cc, (char *)next);

            /*** Normals ***/
            /* prim->n0 is a face normal; prim->n1/2/3 are vertex normals. */
            sprintf(next, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f ",
                prim->n1[0], prim->n1[1], prim->n1[2],
                prim->n2[0], prim->n2[1], prim->n2[2],
                prim->n3[0], prim->n3[1], prim->n3[2]);
            UtilConcatVLA(&normals_str, &norm_str_cc, (char *)next);

            /* Colors */
            /* R, G, B per vertex */
            sprintf(next, "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f ",
                prim->c1[0], prim->c1[1], prim->c1[2],    // vertex 1
                prim->c2[0], prim->c2[1], prim->c2[2],    // vertex 2
                prim->c3[0], prim->c3[1], prim->c3[2]);   // vertex 3
            UtilConcatVLA(&colors_str, &col_str_cc, next);

            /* <p> indices */
            if (TriangleReverse(prim)) {
              sprintf(next, "%i %i %i %i %i %i %i %i %i ",
                  pos, norm, col,
                  pos + 2, norm + 2, col + 2,
                  pos + 1, norm + 1, col + 1);
            } else {
              sprintf(next, "%i %i %i %i %i %i %i %i %i ",
                  pos, norm, col,
                  pos + 1, norm + 1, col + 1,
                  pos + 2, norm + 2, col + 2);
            }
            UtilConcatVLA(&p_str, &p_str_cc, next);

            pos += 3;
            norm += 3;
            col += 3;
            tri += 1;

            free(next);

            if (ReachedXmlNodeSizeLimit(pos_str_cc, norm_str_cc,
                  col_str_cc, p_str_cc)) {

              /* xmlTextWriter uses the LibXML parser module, which enforces a node size
               * limit (10MB as of June 2014). Split oversize nodes into multiple nodes. */
#if COLLADA_DEBUG
              printf("Reached XML_NODE_SIZE_LIMIT (%i) at a=%i, geom=%i\n",
                  XML_NODE_SIZE_LIMIT, a, geom);
#endif

              /* Write the current mesh object and clean up */
              ColladaWriteMeshGeometry(w, geom, pos, positions_str,
                  norm, normals_str, col, colors_str, tri, p_str, geom_mode);

              geom += 1;
              mesh_obj = false;
            }

            break;
          }  // cPrimTriangle

          /* Character and Ellipsoid not implemented */
        case cPrimCharacter:      // 5
          {
#if COLLADA_DEBUG > 1
            printf("Primitive %i: cPrimCharacter\n", a);
#endif
            break;
          }  // cPrimCharacter

        case cPrimEllipsoid:      // 6
          {
#if COLLADA_DEBUG > 1
            printf("Primitive %i: cPrimEllipsoid\n", a);
#endif
            break;
          }  // cPrimEllipsoid

      }  // switch
    }  // for

    xmlTextWriterEndElement(w);  // library_geometries

    VLAFree(positions_str);
    VLAFree(normals_str);
    VLAFree(colors_str);
    VLAFree(p_str);
  }

  /* Effects */
  ColladaWriteLibraryEffects(w, G, trans_len, trans);

  /* Materials */
  ColladaWriteLibraryMaterials(w, trans_len, trans);

  /* Visual Scenes */
  {
    xmlTextWriterStartElement(w, BAD_CAST "library_visual_scenes");

    // TODO: support for multiple scenes (e.g. stored scenes)
    xmlTextWriterStartElement(w, BAD_CAST "visual_scene");
    xmlTextWriterWriteAttribute(w, BAD_CAST "id", BAD_CAST "scene");
    xmlTextWriterWriteAttribute(w, BAD_CAST "name", BAD_CAST "scene");

    /* One node per geometry element. */
    int i, j;
    char *mat = (char *) malloc(100 * sizeof(char));

    for(i = 0; i < geom; i++){
      xmlTextWriterStartElement(w, BAD_CAST "node");
      xmlTextWriterWriteFormatAttribute(w, BAD_CAST "id", "node-geom%i", i);

      /* set center screen pixel as origin */
      float offset[3];
      SceneOriginGet(G, offset);
      invert3f(offset);
      for (j = 0; j < 2; j++) {
        /* camera XY offset (Z handled in camera node) */
        offset[j] += I->Pos[j];
      }
      char *tmp = (char *) malloc(50 * sizeof(char));
      sprintf(tmp, "%5.3f %5.3f %5.3f", offset[0], offset[1], offset[2]);
      xmlTextWriterStartElement(w, BAD_CAST "translate");
      xmlTextWriterWriteAttribute(w, BAD_CAST "sid", BAD_CAST "location");
      xmlTextWriterWriteString(w, BAD_CAST tmp);
      xmlTextWriterEndElement(w);  // translate
      free(tmp);

      xmlTextWriterStartElement(w, BAD_CAST "rotate");
      xmlTextWriterWriteAttribute(w, BAD_CAST "sid", BAD_CAST "rotationZ");
      xmlTextWriterWriteString(w, BAD_CAST "0 0 1 0");
      xmlTextWriterEndElement(w);  // rotate

      xmlTextWriterStartElement(w, BAD_CAST "rotate");
      xmlTextWriterWriteAttribute(w, BAD_CAST "sid", BAD_CAST "rotationY");
      xmlTextWriterWriteString(w, BAD_CAST "0 1 0 0");
      xmlTextWriterEndElement(w);  // rotate

      xmlTextWriterStartElement(w, BAD_CAST "rotate");
      xmlTextWriterWriteAttribute(w, BAD_CAST "sid", BAD_CAST "rotationX");
      xmlTextWriterWriteString(w, BAD_CAST "1 0 0 0");
      xmlTextWriterEndElement(w);  // rotate

      xmlTextWriterStartElement(w, BAD_CAST "scale");
      xmlTextWriterWriteAttribute(w, BAD_CAST "sid", BAD_CAST "scale");
      xmlTextWriterWriteString(w, BAD_CAST "1 1 1");
      xmlTextWriterEndElement(w);  // scale

      xmlTextWriterStartElement(w, BAD_CAST "instance_geometry");
      xmlTextWriterWriteFormatAttribute(w, BAD_CAST "url", "#geom%i", i);
      xmlTextWriterStartElement(w, BAD_CAST "bind_material");
      xmlTextWriterStartElement(w, BAD_CAST "technique_common");
      xmlTextWriterStartElement(w, BAD_CAST "instance_material");
      xmlTextWriterWriteFormatAttribute(w, BAD_CAST "symbol",
          "geom%i-material", i);

      if (TRANS_PRECISION > geom_trans[i]) {
        sprintf(mat, "#default-material");
      }
      else {
        sprintf(mat, "#transparency-%1.2f-material", trans[geom_trans[i]]);
      }
      xmlTextWriterWriteAttribute(w, BAD_CAST "target", BAD_CAST mat);

      xmlTextWriterEndElement(w);  // instance_material
      xmlTextWriterEndElement(w);  // technique_common
      xmlTextWriterEndElement(w);  // bind_material
      xmlTextWriterEndElement(w);  // instance_geometry

      xmlTextWriterEndElement(w);  // node
    }
    free(mat);

    /* background geometry */
    if (!identity && bgbox) {
      xmlTextWriterStartElement(w, BAD_CAST "node");
      xmlTextWriterWriteAttribute(w, BAD_CAST "id", BAD_CAST "node-geom-bg");

      xmlTextWriterStartElement(w, BAD_CAST "translate");
      xmlTextWriterWriteAttribute(w, BAD_CAST "sid", BAD_CAST "location");
      xmlTextWriterWriteString(w, BAD_CAST "0 0 0");
      xmlTextWriterEndElement(w);  // translate

      xmlTextWriterStartElement(w, BAD_CAST "rotate");
      xmlTextWriterWriteAttribute(w, BAD_CAST "sid", BAD_CAST "rotationZ");
      xmlTextWriterWriteString(w, BAD_CAST "0 0 1 0");
      xmlTextWriterEndElement(w);  // rotate

      xmlTextWriterStartElement(w, BAD_CAST "rotate");
      xmlTextWriterWriteAttribute(w, BAD_CAST "sid", BAD_CAST "rotationY");
      xmlTextWriterWriteString(w, BAD_CAST "0 1 0 0");
      xmlTextWriterEndElement(w);  // rotate

      xmlTextWriterStartElement(w, BAD_CAST "rotate");
      xmlTextWriterWriteAttribute(w, BAD_CAST "sid", BAD_CAST "rotationX");
      xmlTextWriterWriteString(w, BAD_CAST "1 0 0 0");
      xmlTextWriterEndElement(w);  // rotate

      xmlTextWriterStartElement(w, BAD_CAST "scale");
      xmlTextWriterWriteAttribute(w, BAD_CAST "sid", BAD_CAST "scale");
      xmlTextWriterWriteString(w, BAD_CAST "1 1 1");
      xmlTextWriterEndElement(w);  // scale

      xmlTextWriterStartElement(w, BAD_CAST "instance_geometry");
      xmlTextWriterWriteAttribute(w, BAD_CAST "url", BAD_CAST "#geom-bg");
      xmlTextWriterStartElement(w, BAD_CAST "bind_material");
      xmlTextWriterStartElement(w, BAD_CAST "technique_common");
      xmlTextWriterStartElement(w, BAD_CAST "instance_material");
      xmlTextWriterWriteAttribute(w, BAD_CAST "symbol",
          BAD_CAST "geom-bg-material");
      xmlTextWriterWriteAttribute(w, BAD_CAST "target",
          BAD_CAST "#bg-material");  // TODO: custom materials
      xmlTextWriterEndElement(w);  // instance_material
      xmlTextWriterEndElement(w);  // technique_common
      xmlTextWriterEndElement(w);  // bind_material
      xmlTextWriterEndElement(w);  // instance_geometry

      xmlTextWriterEndElement(w);  // node
    }

    /* Camera node */
    char *matrix = (char *) malloc(400 * sizeof(char));  // enough for 16 floats
    if (!identity) {
      float m1[16];

      // m1 = RotMatrix * T(0, 0, -Pos[z])
      identity44f(m1);

      // camera z
      m1[11] = -I->Pos[2];

      // RotMatrix
      left_multiply44f44f(SceneGetMatrix(G), m1);

      sprintf(matrix,
          "\n%5.3f %5.3f %5.3f %5.3f "
          "\n%5.3f %5.3f %5.3f %5.3f "
          "\n%5.3f %5.3f %5.3f %5.3f "
          "\n%5.3f %5.3f %5.3f %5.3f ",
          m1[0], m1[1], m1[2], m1[3],
          m1[4], m1[5], m1[6], m1[7],
          m1[8], m1[9], m1[10], m1[11],
          m1[12], m1[13], m1[14], m1[15]);

      xmlTextWriterStartElement(w, BAD_CAST "node");
      xmlTextWriterWriteAttribute(w, BAD_CAST "id", BAD_CAST "node-camera");

      xmlTextWriterWriteElement(w, BAD_CAST "matrix", BAD_CAST matrix);

      xmlTextWriterStartElement(w, BAD_CAST "scale");
      xmlTextWriterWriteAttribute(w, BAD_CAST "sid", BAD_CAST "scale");
      xmlTextWriterWriteString(w, BAD_CAST "1 1 1");
      xmlTextWriterEndElement(w);  // scale

      xmlTextWriterStartElement(w, BAD_CAST "instance_camera");
      xmlTextWriterWriteAttribute(w, BAD_CAST "url", BAD_CAST "#camera");
      xmlTextWriterEndElement(w);  // instance_camera

      xmlTextWriterEndElement(w);  // node

    }

    /* Light nodes */
    if (!identity && lighting) {
      if(lc > 0){

        /* First light is the ambient light.  For light_count of 2-10,
         * we add a directional light, e.g. light, light2, etc. */

        int i;

        xmlTextWriterStartElement(w, BAD_CAST "node");
        xmlTextWriterWriteFormatAttribute(w, BAD_CAST "id", "node-ambient-light");
        xmlTextWriterWriteAttribute(w, BAD_CAST "type", BAD_CAST "NODE");
        xmlTextWriterStartElement(w, BAD_CAST "instance_light");
        xmlTextWriterWriteAttribute(w, BAD_CAST "url", BAD_CAST "#ambient-light");
        xmlTextWriterEndElement(w);  // instance_light

        xmlTextWriterEndElement(w);  // node

        const float *light_pos_ptr;
        float light_pos[3];
        for(i = 1; i < lc; i++){
          switch(i){
            case 1:
              light_pos_ptr = SettingGetGlobal_3fv(G, cSetting_light);
              break;
            case 2:
              light_pos_ptr = SettingGetGlobal_3fv(G, cSetting_light2);
              break;
            case 3:
              light_pos_ptr = SettingGetGlobal_3fv(G, cSetting_light3);
              break;
            case 4:
              light_pos_ptr = SettingGetGlobal_3fv(G, cSetting_light4);
              break;
            case 5:
              light_pos_ptr = SettingGetGlobal_3fv(G, cSetting_light5);
              break;
            case 6:
              light_pos_ptr = SettingGetGlobal_3fv(G, cSetting_light6);
              break;
            case 7:
              light_pos_ptr = SettingGetGlobal_3fv(G, cSetting_light7);
              break;
            case 8:
              light_pos_ptr = SettingGetGlobal_3fv(G, cSetting_light8);
              break;
            default:
              light_pos_ptr = SettingGetGlobal_3fv(G, cSetting_light9);
              break;
          }

          copy3f(light_pos_ptr, light_pos);
          normalize3f(light_pos);

          xmlTextWriterStartElement(w, BAD_CAST "node");
          xmlTextWriterWriteFormatAttribute(w, BAD_CAST "id", "node-pymol-light%i", i);
          xmlTextWriterWriteAttribute(w, BAD_CAST "type", BAD_CAST "NODE");

          char *lookat = (char *) malloc(sizeof(char) * 50);
          sprintf(lookat,
              "%5.3f %5.3f %5.3f "  // position of light on unit sphere
              "0 0 0 "  // pointed toward origin
              "0 1 0",  // positive y axis up (arbitrary)
              /* negative because PyMOL gives direction of light, not position */
              -light_pos[0], -light_pos[1], -light_pos[2]);

          xmlTextWriterWriteElement(w, BAD_CAST "lookat", BAD_CAST lookat);

          xmlTextWriterWriteElement(w, BAD_CAST "matrix", BAD_CAST matrix);

          xmlTextWriterStartElement(w, BAD_CAST "instance_light");
          xmlTextWriterWriteAttribute(w, BAD_CAST "url", BAD_CAST "#pymol-light");
          xmlTextWriterEndElement(w);  // instance_light

          xmlTextWriterEndElement(w);  // node
        }
      }
    }

    free(matrix);

    xmlTextWriterEndElement(w);  // visual_scene
    xmlTextWriterEndElement(w);  // library_visual_scenes
  }

  /* Scene */
  {
    xmlTextWriterStartElement(w, BAD_CAST "scene");
    xmlTextWriterStartElement(w, BAD_CAST "instance_visual_scene");
    xmlTextWriterWriteAttribute(w, BAD_CAST "url", BAD_CAST "#scene");
    xmlTextWriterEndElement(w);  // instance_visual_scene
    xmlTextWriterEndElement(w);  // scene
  }


  xmlTextWriterEndElement(w);  // COLLADA

  /* Close the document */
  rc = xmlTextWriterEndDocument(w);
  if (rc < 0) {
    printf("ColladaRender: Error at xmlTextWriterEndDocument.\n");
    // return;
  }

  /*
   * Dump the document to buffer
   */
  xmlDocDumpFormatMemory(doc, &tmp, &buffersize, 1);

  /* Concat buffer to VLA output */
  UtilConcatVLA(&vla, &cc, (char *) tmp);

  /*
   * Free associated memory.
   */
  free(trans);
  free(geom_trans);

  xmlFree(tmp);
  xmlFreeTextWriter(w);
  xmlFreeDoc(doc);

  *vla_ptr = vla;

#endif // _HAVE_LIBXML
}
