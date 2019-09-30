/*
OpenVR for PyMOL Copyright Notice
=====================================

The OpenVR for PyMOL source code is copyrighted, but you can freely use and
copy it as long as you don't change or remove any of the Copyright notices.
OpenVR for PyMOL is made available under the following open-source license
terms:

------------------------------------------------------------------------------
Copyright (c) 2018 EPAM Systems, Inc.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
------------------------------------------------------------------------------

*/

#ifndef _H_OpenVRStub
#define _H_OpenVRStub

// system headers
#include "openvr.h"

namespace vr {
namespace stub {

void VR_StubEnable(bool on = true);
bool VR_IsStubEnabled(void);

/** Finds the active installation of the VR API and initializes it. The provided path must be absolute
* or relative to the current working directory. These are the local install versions of the equivalent
* functions in steamvr.h and will work without a local Steam install.
*
* This path is to the "root" of the VR API install. That's the directory with
* the "drivers" directory and a platform (i.e. "win32") directory in it, not the directory with the DLL itself.
*
* pStartupInfo is reserved for future use.
*/
IVRSystem *VR_Init(EVRInitError *peError, EVRApplicationType eApplicationType, const char *pStartupInfo = nullptr);

/** unloads vrclient.dll. Any interface pointers from the interface are
* invalid after this point */
void VR_Shutdown();

/** Returns true if there is an HMD attached. This check is as lightweight as possible and
* can be called outside of VR_Init/VR_Shutdown. It should be used when an application wants
* to know if initializing VR is a possibility but isn't ready to take that step yet.
*/
bool VR_IsHmdPresent();

/** Returns an English string for an EVRInitError. Applications should call VR_GetVRInitErrorAsSymbol instead and
* use that as a key to look up their own localized error message. This function may be called outside of VR_Init()/VR_Shutdown(). */
const char *VR_GetVRInitErrorAsEnglishDescription(EVRInitError error);

IVRSystem *VRSystem();
IVRCompositor *VRCompositor();
IVRInput *VRInput();

#define VRSTUB_NOT_IMPLEMENTED \
  throw __FUNCTION__; // " not implemented";

class VRSystemStub: public IVRSystem {
public:
  // ------------------------------------
  // Display Methods
  // ------------------------------------

  /** Suggested size for the intermediate render target that the distortion pulls from. */
  virtual void GetRecommendedRenderTargetSize( uint32_t *pnWidth, uint32_t *pnHeight );

  /** The projection matrix for the specified eye */
  virtual HmdMatrix44_t GetProjectionMatrix( EVREye eEye, float fNearZ, float fFarZ ) { VRSTUB_NOT_IMPLEMENTED }

  /** The components necessary to build your own projection matrix in case your
  * application is doing something fancy like infinite Z */
  virtual void GetProjectionRaw( EVREye eEye, float *pfLeft, float *pfRight, float *pfTop, float *pfBottom );

  /** Gets the result of the distortion function for the specified eye and input UVs. UVs go from 0,0 in 
  * the upper left of that eye's viewport and 1,1 in the lower right of that eye's viewport.
  * Returns true for success. Otherwise, returns false, and distortion coordinates are not suitable. */
  virtual bool ComputeDistortion( EVREye eEye, float fU, float fV, DistortionCoordinates_t *pDistortionCoordinates ) { VRSTUB_NOT_IMPLEMENTED }

  /** Returns the transform from eye space to the head space. Eye space is the per-eye flavor of head
  * space that provides stereo disparity. Instead of Model * View * Projection the sequence is Model * View * Eye^-1 * Projection. 
  * Normally View and Eye^-1 will be multiplied together and treated as View in your application. 
  */
  virtual HmdMatrix34_t GetEyeToHeadTransform( EVREye eEye );

  /** Returns the number of elapsed seconds since the last recorded vsync event. This 
  *	will come from a vsync timer event in the timer if possible or from the application-reported
  *   time if that is not available. If no vsync times are available the function will 
  *   return zero for vsync time and frame counter and return false from the method. */
  virtual bool GetTimeSinceLastVsync( float *pfSecondsSinceLastVsync, uint64_t *pulFrameCounter ) { VRSTUB_NOT_IMPLEMENTED }

  /** [D3D9 Only]
  * Returns the adapter index that the user should pass into CreateDevice to set up D3D9 in such
  * a way that it can go full screen exclusive on the HMD. Returns -1 if there was an error.
  */
  virtual int32_t GetD3D9AdapterIndex() { VRSTUB_NOT_IMPLEMENTED }

  /** [D3D10/11 Only]
  * Returns the adapter index that the user should pass into EnumAdapters to create the device 
  * and swap chain in DX10 and DX11. If an error occurs the index will be set to -1.
  */
  virtual void GetDXGIOutputInfo( int32_t *pnAdapterIndex ) { VRSTUB_NOT_IMPLEMENTED }

  /**
   * Returns platform- and texture-type specific adapter identification so that applications and the
   * compositor are creating textures and swap chains on the same GPU. If an error occurs the device
   * will be set to 0.
   * pInstance is an optional parameter that is required only when textureType is TextureType_Vulkan.
   * [D3D10/11/12 Only (D3D9 Not Supported)]
   *  Returns the adapter LUID that identifies the GPU attached to the HMD. The user should
   *  enumerate all adapters using IDXGIFactory::EnumAdapters and IDXGIAdapter::GetDesc to find
   *  the adapter with the matching LUID, or use IDXGIFactory4::EnumAdapterByLuid.
   *  The discovered IDXGIAdapter should be used to create the device and swap chain.
   * [Vulkan Only]
   *  Returns the VkPhysicalDevice that should be used by the application.
   *  pInstance must be the instance the application will use to query for the VkPhysicalDevice.  The application
   *  must create the VkInstance with extensions returned by IVRCompositor::GetVulkanInstanceExtensionsRequired enabled.
   * [macOS Only]
   *  For TextureType_IOSurface returns the id<MTLDevice> that should be used by the application.
   *  On 10.13+ for TextureType_OpenGL returns the 'registryId' of the renderer which should be used
   *   by the application. See Apple Technical Q&A QA1168 for information on enumerating GL Renderers, and the
   *   new kCGLRPRegistryIDLow and kCGLRPRegistryIDHigh CGLRendererProperty values in the 10.13 SDK.
   *  Pre 10.13 for TextureType_OpenGL returns 0, as there is no dependable way to correlate the HMDs MTLDevice
   *   with a GL Renderer.
   */
  virtual void GetOutputDevice( uint64_t *pnDevice, ETextureType textureType, VkInstance_T *pInstance = nullptr ) { VRSTUB_NOT_IMPLEMENTED }

  // ------------------------------------
  // Display Mode methods
  // ------------------------------------

  /** Use to determine if the headset display is part of the desktop (i.e. extended) or hidden (i.e. direct mode). */
  virtual bool IsDisplayOnDesktop() { VRSTUB_NOT_IMPLEMENTED }

  /** Set the display visibility (true = extended, false = direct mode).  Return value of true indicates that the change was successful. */
  virtual bool SetDisplayVisibility( bool bIsVisibleOnDesktop ) { VRSTUB_NOT_IMPLEMENTED }

  // ------------------------------------
  // Tracking Methods
  // ------------------------------------

  /** The pose that the tracker thinks that the HMD will be in at the specified number of seconds into the 
  * future. Pass 0 to get the state at the instant the method is called. Most of the time the application should
  * calculate the time until the photons will be emitted from the display and pass that time into the method.
  *
  * This is roughly analogous to the inverse of the view matrix in most applications, though 
  * many games will need to do some additional rotation or translation on top of the rotation
  * and translation provided by the head pose.
  *
  * For devices where bPoseIsValid is true the application can use the pose to position the device
  * in question. The provided array can be any size up to k_unMaxTrackedDeviceCount. 
  *
  * Seated experiences should call this method with TrackingUniverseSeated and receive poses relative
  * to the seated zero pose. Standing experiences should call this method with TrackingUniverseStanding 
  * and receive poses relative to the Chaperone Play Area. TrackingUniverseRawAndUncalibrated should 
  * probably not be used unless the application is the Chaperone calibration tool itself, but will provide
  * poses relative to the hardware-specific coordinate system in the driver.
  */
  virtual void GetDeviceToAbsoluteTrackingPose( ETrackingUniverseOrigin eOrigin, float fPredictedSecondsToPhotonsFromNow, VR_ARRAY_COUNT(unTrackedDevicePoseArrayCount) TrackedDevicePose_t *pTrackedDevicePoseArray, uint32_t unTrackedDevicePoseArrayCount ) { VRSTUB_NOT_IMPLEMENTED }

  /** Sets the zero pose for the seated tracker coordinate system to the current position and yaw of the HMD. After 
  * ResetSeatedZeroPose all GetDeviceToAbsoluteTrackingPose calls that pass TrackingUniverseSeated as the origin 
  * will be relative to this new zero pose. The new zero coordinate system will not change the fact that the Y axis 
  * is up in the real world, so the next pose returned from GetDeviceToAbsoluteTrackingPose after a call to 
  * ResetSeatedZeroPose may not be exactly an identity matrix.
  *
  * NOTE: This function overrides the user's previously saved seated zero pose and should only be called as the result of a user action. 
  * Users are also able to set their seated zero pose via the OpenVR Dashboard.
  **/
  virtual void ResetSeatedZeroPose() { VRSTUB_NOT_IMPLEMENTED }

  /** Returns the transform from the seated zero pose to the standing absolute tracking system. This allows 
  * applications to represent the seated origin to used or transform object positions from one coordinate
  * system to the other. 
  *
  * The seated origin may or may not be inside the Play Area or Collision Bounds returned by IVRChaperone. Its position 
  * depends on what the user has set from the Dashboard settings and previous calls to ResetSeatedZeroPose. */
  virtual HmdMatrix34_t GetSeatedZeroPoseToStandingAbsoluteTrackingPose() { VRSTUB_NOT_IMPLEMENTED }

  /** Returns the transform from the tracking origin to the standing absolute tracking system. This allows
  * applications to convert from raw tracking space to the calibrated standing coordinate system. */
  virtual HmdMatrix34_t GetRawZeroPoseToStandingAbsoluteTrackingPose() { VRSTUB_NOT_IMPLEMENTED }

  /** Get a sorted array of device indices of a given class of tracked devices (e.g. controllers).  Devices are sorted right to left
  * relative to the specified tracked device (default: hmd -- pass in -1 for absolute tracking space).  Returns the number of devices
  * in the list, or the size of the array needed if not large enough. */
  virtual uint32_t GetSortedTrackedDeviceIndicesOfClass( ETrackedDeviceClass eTrackedDeviceClass, VR_ARRAY_COUNT(unTrackedDeviceIndexArrayCount) vr::TrackedDeviceIndex_t *punTrackedDeviceIndexArray, uint32_t unTrackedDeviceIndexArrayCount, vr::TrackedDeviceIndex_t unRelativeToTrackedDeviceIndex = k_unTrackedDeviceIndex_Hmd ) { VRSTUB_NOT_IMPLEMENTED }

  /** Returns the level of activity on the device. */
  virtual EDeviceActivityLevel GetTrackedDeviceActivityLevel( vr::TrackedDeviceIndex_t unDeviceId ) { VRSTUB_NOT_IMPLEMENTED }

  /** Convenience utility to apply the specified transform to the specified pose.
  *   This properly transforms all pose components, including velocity and angular velocity
  */
  virtual void ApplyTransform( TrackedDevicePose_t *pOutputPose, const TrackedDevicePose_t *pTrackedDevicePose, const HmdMatrix34_t *pTransform ) { VRSTUB_NOT_IMPLEMENTED }

  /** Returns the device index associated with a specific role, for example the left hand or the right hand. */
  virtual vr::TrackedDeviceIndex_t GetTrackedDeviceIndexForControllerRole( vr::ETrackedControllerRole unDeviceType ) { VRSTUB_NOT_IMPLEMENTED }

  /** Returns the controller type associated with a device index. */
  virtual vr::ETrackedControllerRole GetControllerRoleForTrackedDeviceIndex( vr::TrackedDeviceIndex_t unDeviceIndex ) { VRSTUB_NOT_IMPLEMENTED }

  // ------------------------------------
  // Property methods
  // ------------------------------------

  /** Returns the device class of a tracked device. If there has not been a device connected in this slot
  * since the application started this function will return TrackedDevice_Invalid. For previous detected
  * devices the function will return the previously observed device class. 
  *
  * To determine which devices exist on the system, just loop from 0 to k_unMaxTrackedDeviceCount and check
  * the device class. Every device with something other than TrackedDevice_Invalid is associated with an 
  * actual tracked device. */
  virtual ETrackedDeviceClass GetTrackedDeviceClass( vr::TrackedDeviceIndex_t unDeviceIndex );

  /** Returns true if there is a device connected in this slot. */
  virtual bool IsTrackedDeviceConnected( vr::TrackedDeviceIndex_t unDeviceIndex ) { VRSTUB_NOT_IMPLEMENTED }

  /** Returns a bool property. If the device index is not valid or the property is not a bool type this function will return false. */
  virtual bool GetBoolTrackedDeviceProperty( vr::TrackedDeviceIndex_t unDeviceIndex, ETrackedDeviceProperty prop, ETrackedPropertyError *pError = 0L ) { VRSTUB_NOT_IMPLEMENTED }

  /** Returns a float property. If the device index is not valid or the property is not a float type this function will return 0. */
  virtual float GetFloatTrackedDeviceProperty( vr::TrackedDeviceIndex_t unDeviceIndex, ETrackedDeviceProperty prop, ETrackedPropertyError *pError = 0L ) { VRSTUB_NOT_IMPLEMENTED }

  /** Returns an int property. If the device index is not valid or the property is not a int type this function will return 0. */
  virtual int32_t GetInt32TrackedDeviceProperty( vr::TrackedDeviceIndex_t unDeviceIndex, ETrackedDeviceProperty prop, ETrackedPropertyError *pError = 0L ) { VRSTUB_NOT_IMPLEMENTED }

  /** Returns a uint64 property. If the device index is not valid or the property is not a uint64 type this function will return 0. */
  virtual uint64_t GetUint64TrackedDeviceProperty( vr::TrackedDeviceIndex_t unDeviceIndex, ETrackedDeviceProperty prop, ETrackedPropertyError *pError = 0L ) { VRSTUB_NOT_IMPLEMENTED }

  /** Returns a matrix property. If the device index is not valid or the property is not a matrix type, this function will return identity. */
  virtual HmdMatrix34_t GetMatrix34TrackedDeviceProperty( vr::TrackedDeviceIndex_t unDeviceIndex, ETrackedDeviceProperty prop, ETrackedPropertyError *pError = 0L ) { VRSTUB_NOT_IMPLEMENTED }

  /** Returns an array of one type of property. If the device index is not valid or the property is not a single value or an array of the specified type,
  * this function will return 0. Otherwise it returns the number of bytes necessary to hold the array of properties. If unBufferSize is
  * greater than the returned size and pBuffer is non-NULL, pBuffer is filled with the contents of array of properties. */
  virtual uint32_t GetArrayTrackedDeviceProperty( vr::TrackedDeviceIndex_t unDeviceIndex, ETrackedDeviceProperty prop, PropertyTypeTag_t propType, void *pBuffer, uint32_t unBufferSize, ETrackedPropertyError *pError = 0L ) { VRSTUB_NOT_IMPLEMENTED }

  /** Returns a string property. If the device index is not valid or the property is not a string type this function will 
  * return 0. Otherwise it returns the length of the number of bytes necessary to hold this string including the trailing
  * null. Strings will always fit in buffers of k_unMaxPropertyStringSize characters. */
  virtual uint32_t GetStringTrackedDeviceProperty( vr::TrackedDeviceIndex_t unDeviceIndex, ETrackedDeviceProperty prop, VR_OUT_STRING() char *pchValue, uint32_t unBufferSize, ETrackedPropertyError *pError = 0L );

  /** returns a string that corresponds with the specified property error. The string will be the name 
  * of the error enum value for all valid error codes */
  virtual const char *GetPropErrorNameFromEnum( ETrackedPropertyError error ) { VRSTUB_NOT_IMPLEMENTED }

  // ------------------------------------
  // Event methods
  // ------------------------------------

  /** Returns true and fills the event with the next event on the queue if there is one. If there are no events
  * this method returns false. uncbVREvent should be the size in bytes of the VREvent_t struct */
  virtual bool PollNextEvent( VREvent_t *pEvent, uint32_t uncbVREvent );

  /** Returns true and fills the event with the next event on the queue if there is one. If there are no events
  * this method returns false. Fills in the pose of the associated tracked device in the provided pose struct. 
  * This pose will always be older than the call to this function and should not be used to render the device. 
  uncbVREvent should be the size in bytes of the VREvent_t struct */
  virtual bool PollNextEventWithPose( ETrackingUniverseOrigin eOrigin, VREvent_t *pEvent, uint32_t uncbVREvent, vr::TrackedDevicePose_t *pTrackedDevicePose ) { VRSTUB_NOT_IMPLEMENTED }

  /** returns the name of an EVREvent enum value */
  virtual const char *GetEventTypeNameFromEnum( EVREventType eType ) { VRSTUB_NOT_IMPLEMENTED }

  // ------------------------------------
  // Rendering helper methods
  // ------------------------------------

  /** Returns the hidden area mesh for the current HMD. The pixels covered by this mesh will never be seen by the user after the lens distortion is
  * applied based on visibility to the panels. If this HMD does not have a hidden area mesh, the vertex data and count will be NULL and 0 respectively.
  * This mesh is meant to be rendered into the stencil buffer (or into the depth buffer setting nearz) before rendering each eye's view. 
  * This will improve performance by letting the GPU early-reject pixels the user will never see before running the pixel shader.
  * NOTE: Render this mesh with backface culling disabled since the winding order of the vertices can be different per-HMD or per-eye.
  * Setting the bInverse argument to true will produce the visible area mesh that is commonly used in place of full-screen quads. The visible area mesh covers all of the pixels the hidden area mesh does not cover.
  * Setting the bLineLoop argument will return a line loop of vertices in HiddenAreaMesh_t->pVertexData with HiddenAreaMesh_t->unTriangleCount set to the number of vertices.
  */
  virtual HiddenAreaMesh_t GetHiddenAreaMesh( EVREye eEye, EHiddenAreaMeshType type = k_eHiddenAreaMesh_Standard ) { VRSTUB_NOT_IMPLEMENTED }

  // ------------------------------------
  // Controller methods
  // ------------------------------------

  /** Fills the supplied struct with the current state of the controller. Returns false if the controller index
  * is invalid. */
  virtual bool GetControllerState( vr::TrackedDeviceIndex_t unControllerDeviceIndex, vr::VRControllerState_t *pControllerState, uint32_t unControllerStateSize ) { VRSTUB_NOT_IMPLEMENTED }

  /** fills the supplied struct with the current state of the controller and the provided pose with the pose of 
  * the controller when the controller state was updated most recently. Use this form if you need a precise controller
  * pose as input to your application when the user presses or releases a button. */
  virtual bool GetControllerStateWithPose( ETrackingUniverseOrigin eOrigin, vr::TrackedDeviceIndex_t unControllerDeviceIndex, vr::VRControllerState_t *pControllerState, uint32_t unControllerStateSize, TrackedDevicePose_t *pTrackedDevicePose ) { VRSTUB_NOT_IMPLEMENTED }

  /** Trigger a single haptic pulse on a controller. After this call the application may not trigger another haptic pulse on this controller
  * and axis combination for 5ms. */
  virtual void TriggerHapticPulse( vr::TrackedDeviceIndex_t unControllerDeviceIndex, uint32_t unAxisId, unsigned short usDurationMicroSec ) { VRSTUB_NOT_IMPLEMENTED }

  /** returns the name of an EVRButtonId enum value */
  virtual const char *GetButtonIdNameFromEnum( EVRButtonId eButtonId ) { VRSTUB_NOT_IMPLEMENTED }

  /** returns the name of an EVRControllerAxisType enum value */
  virtual const char *GetControllerAxisTypeNameFromEnum( EVRControllerAxisType eAxisType ) { VRSTUB_NOT_IMPLEMENTED }

  /** Returns true if this application is receiving input from the system. This would return false if 
  * system-related functionality is consuming the input stream. */
  virtual bool IsInputAvailable() { VRSTUB_NOT_IMPLEMENTED }

  /** Returns true SteamVR is drawing controllers on top of the application. Applications should consider
  * not drawing anything attached to the user's hands in this case. */
  virtual bool IsSteamVRDrawingControllers() { VRSTUB_NOT_IMPLEMENTED }

  /** Returns true if the user has put SteamVR into a mode that is distracting them from the application.
  * For applications where this is appropriate, the application should pause ongoing activity. */
  virtual bool ShouldApplicationPause() { VRSTUB_NOT_IMPLEMENTED }

  /** Returns true if SteamVR is doing significant rendering work and the game should do what it can to reduce
  * its own workload. One common way to do this is to reduce the size of the render target provided for each eye. */
  virtual bool ShouldApplicationReduceRenderingWork() { VRSTUB_NOT_IMPLEMENTED }

  // ------------------------------------
  // Debug Methods
  // ------------------------------------

  /** Sends a request to the driver for the specified device and returns the response. The maximum response size is 32k,
  * but this method can be called with a smaller buffer. If the response exceeds the size of the buffer, it is truncated. 
  * The size of the response including its terminating null is returned. */
  virtual uint32_t DriverDebugRequest( vr::TrackedDeviceIndex_t unDeviceIndex, const char *pchRequest, VR_OUT_STRING() char *pchResponseBuffer, uint32_t unResponseBufferSize ) { VRSTUB_NOT_IMPLEMENTED }

  // ------------------------------------
  // Firmware methods
  // ------------------------------------

  /** Performs the actual firmware update if applicable. 
   * The following events will be sent, if VRFirmwareError_None was returned: VREvent_FirmwareUpdateStarted, VREvent_FirmwareUpdateFinished 
   * Use the properties Prop_Firmware_UpdateAvailable_Bool, Prop_Firmware_ManualUpdate_Bool, and Prop_Firmware_ManualUpdateURL_String
   * to figure our whether a firmware update is available, and to figure out whether its a manual update 
   * Prop_Firmware_ManualUpdateURL_String should point to an URL describing the manual update process */
  virtual vr::EVRFirmwareError PerformFirmwareUpdate( vr::TrackedDeviceIndex_t unDeviceIndex ) { VRSTUB_NOT_IMPLEMENTED }

  // ------------------------------------
  // Application life cycle methods
  // ------------------------------------

  /** Call this to acknowledge to the system that VREvent_Quit has been received and that the process is exiting.
  * This extends the timeout until the process is killed. */
  virtual void AcknowledgeQuit_Exiting() { VRSTUB_NOT_IMPLEMENTED }

  /** Call this to tell the system that the user is being prompted to save data. This
  * halts the timeout and dismisses the dashboard (if it was up). Applications should be sure to actually 
  * prompt the user to save and then exit afterward, otherwise the user will be left in a confusing state. */
  virtual void AcknowledgeQuit_UserPrompt() { VRSTUB_NOT_IMPLEMENTED }

};

class VRCompositorStub: public IVRCompositor {
public:
  /** Sets tracking space returned by WaitGetPoses */
  virtual void SetTrackingSpace( ETrackingUniverseOrigin eOrigin ) { VRSTUB_NOT_IMPLEMENTED }

  /** Gets current tracking space returned by WaitGetPoses */
  virtual ETrackingUniverseOrigin GetTrackingSpace() { VRSTUB_NOT_IMPLEMENTED }

  /** Scene applications should call this function to get poses to render with (and optionally poses predicted an additional frame out to use for gameplay).
  * This function will block until "running start" milliseconds before the start of the frame, and should be called at the last moment before needing to
  * start rendering.
  *
  * Return codes:
  *	- IsNotSceneApplication (make sure to call VR_Init with VRApplicaiton_Scene)
  *	- DoNotHaveFocus (some other app has taken focus - this will throttle the call to 10hz to reduce the impact on that app)
  */
  virtual EVRCompositorError WaitGetPoses( VR_ARRAY_COUNT(unRenderPoseArrayCount) TrackedDevicePose_t* pRenderPoseArray, uint32_t unRenderPoseArrayCount,
    VR_ARRAY_COUNT(unGamePoseArrayCount) TrackedDevicePose_t* pGamePoseArray, uint32_t unGamePoseArrayCount );

  /** Get the last set of poses returned by WaitGetPoses. */
  virtual EVRCompositorError GetLastPoses( VR_ARRAY_COUNT( unRenderPoseArrayCount ) TrackedDevicePose_t* pRenderPoseArray, uint32_t unRenderPoseArrayCount,
    VR_ARRAY_COUNT( unGamePoseArrayCount ) TrackedDevicePose_t* pGamePoseArray, uint32_t unGamePoseArrayCount ) { VRSTUB_NOT_IMPLEMENTED }

  /** Interface for accessing last set of poses returned by WaitGetPoses one at a time.
  * Returns VRCompositorError_IndexOutOfRange if unDeviceIndex not less than k_unMaxTrackedDeviceCount otherwise VRCompositorError_None.
  * It is okay to pass NULL for either pose if you only want one of the values. */
  virtual EVRCompositorError GetLastPoseForTrackedDeviceIndex( TrackedDeviceIndex_t unDeviceIndex, TrackedDevicePose_t *pOutputPose, TrackedDevicePose_t *pOutputGamePose ) { VRSTUB_NOT_IMPLEMENTED }

  /** Updated scene texture to display. If pBounds is NULL the entire texture will be used.  If called from an OpenGL app, consider adding a glFlush after
  * Submitting both frames to signal the driver to start processing, otherwise it may wait until the command buffer fills up, causing the app to miss frames.
  *
  * OpenGL dirty state:
  *	glBindTexture
  *
  * Return codes:
  *	- IsNotSceneApplication (make sure to call VR_Init with VRApplicaiton_Scene)
  *	- DoNotHaveFocus (some other app has taken focus)
  *	- TextureIsOnWrongDevice (application did not use proper AdapterIndex - see IVRSystem.GetDXGIOutputInfo)
  *	- SharedTexturesNotSupported (application needs to call CreateDXGIFactory1 or later before creating DX device)
  *	- TextureUsesUnsupportedFormat (scene textures must be compatible with DXGI sharing rules - e.g. uncompressed, no mips, etc.)
  *	- InvalidTexture (usually means bad arguments passed in)
  *	- AlreadySubmitted (app has submitted two left textures or two right textures in a single frame - i.e. before calling WaitGetPoses again)
  */
  virtual EVRCompositorError Submit( EVREye eEye, const Texture_t *pTexture, const VRTextureBounds_t* pBounds = 0, EVRSubmitFlags nSubmitFlags = Submit_Default );

  /** Clears the frame that was sent with the last call to Submit. This will cause the 
  * compositor to show the grid until Submit is called again. */
  virtual void ClearLastSubmittedFrame() { VRSTUB_NOT_IMPLEMENTED }

  /** Call immediately after presenting your app's window (i.e. companion window) to unblock the compositor.
  * This is an optional call, which only needs to be used if you can't instead call WaitGetPoses immediately after Present.
  * For example, if your engine's render and game loop are not on separate threads, or blocking the render thread until 3ms before the next vsync would
  * introduce a deadlock of some sort.  This function tells the compositor that you have finished all rendering after having Submitted buffers for both
  * eyes, and it is free to start its rendering work.  This should only be called from the same thread you are rendering on. */
  virtual void PostPresentHandoff() { VRSTUB_NOT_IMPLEMENTED }

  /** Returns true if timing data is filled it.  Sets oldest timing info if nFramesAgo is larger than the stored history.
  * Be sure to set timing.size = sizeof(Compositor_FrameTiming) on struct passed in before calling this function. */
  virtual bool GetFrameTiming( Compositor_FrameTiming *pTiming, uint32_t unFramesAgo = 0 ) { VRSTUB_NOT_IMPLEMENTED }

  /** Interface for copying a range of timing data.  Frames are returned in ascending order (oldest to newest) with the last being the most recent frame.
  * Only the first entry's m_nSize needs to be set, as the rest will be inferred from that.  Returns total number of entries filled out. */
  virtual uint32_t GetFrameTimings( Compositor_FrameTiming *pTiming, uint32_t nFrames ) { VRSTUB_NOT_IMPLEMENTED }

  /** Returns the time in seconds left in the current (as identified by FrameTiming's frameIndex) frame.
  * Due to "running start", this value may roll over to the next frame before ever reaching 0.0. */
  virtual float GetFrameTimeRemaining() { VRSTUB_NOT_IMPLEMENTED }

  /** Fills out stats accumulated for the last connected application.  Pass in sizeof( Compositor_CumulativeStats ) as second parameter. */
  virtual void GetCumulativeStats( Compositor_CumulativeStats *pStats, uint32_t nStatsSizeInBytes ) { VRSTUB_NOT_IMPLEMENTED }

  /** Fades the view on the HMD to the specified color. The fade will take fSeconds, and the color values are between
  * 0.0 and 1.0. This color is faded on top of the scene based on the alpha parameter. Removing the fade color instantly 
  * would be FadeToColor( 0.0, 0.0, 0.0, 0.0, 0.0 ).  Values are in un-premultiplied alpha space. */
  virtual void FadeToColor( float fSeconds, float fRed, float fGreen, float fBlue, float fAlpha, bool bBackground = false ) { VRSTUB_NOT_IMPLEMENTED }

  /** Get current fade color value. */
  virtual HmdColor_t GetCurrentFadeColor( bool bBackground = false ) { VRSTUB_NOT_IMPLEMENTED }

  /** Fading the Grid in or out in fSeconds */
  virtual void FadeGrid( float fSeconds, bool bFadeIn ) { VRSTUB_NOT_IMPLEMENTED }

  /** Get current alpha value of grid. */
  virtual float GetCurrentGridAlpha() { VRSTUB_NOT_IMPLEMENTED }

  /** Override the skybox used in the compositor (e.g. for during level loads when the app can't feed scene images fast enough)
  * Order is Front, Back, Left, Right, Top, Bottom.  If only a single texture is passed, it is assumed in lat-long format.
  * If two are passed, it is assumed a lat-long stereo pair. */
  virtual EVRCompositorError SetSkyboxOverride( VR_ARRAY_COUNT( unTextureCount ) const Texture_t *pTextures, uint32_t unTextureCount ) { VRSTUB_NOT_IMPLEMENTED }

  /** Resets compositor skybox back to defaults. */
  virtual void ClearSkyboxOverride() { VRSTUB_NOT_IMPLEMENTED }

  /** Brings the compositor window to the front. This is useful for covering any other window that may be on the HMD
  * and is obscuring the compositor window. */
  virtual void CompositorBringToFront() { VRSTUB_NOT_IMPLEMENTED }

  /** Pushes the compositor window to the back. This is useful for allowing other applications to draw directly to the HMD. */
  virtual void CompositorGoToBack() { VRSTUB_NOT_IMPLEMENTED }

  /** Tells the compositor process to clean up and exit. You do not need to call this function at shutdown. Under normal 
  * circumstances the compositor will manage its own life cycle based on what applications are running. */
  virtual void CompositorQuit() { VRSTUB_NOT_IMPLEMENTED }

  /** Return whether the compositor is fullscreen */
  virtual bool IsFullscreen() { VRSTUB_NOT_IMPLEMENTED }

  /** Returns the process ID of the process that is currently rendering the scene */
  virtual uint32_t GetCurrentSceneFocusProcess() { VRSTUB_NOT_IMPLEMENTED }

  /** Returns the process ID of the process that rendered the last frame (or 0 if the compositor itself rendered the frame.)
  * Returns 0 when fading out from an app and the app's process Id when fading into an app. */
  virtual uint32_t GetLastFrameRenderer() { VRSTUB_NOT_IMPLEMENTED }

  /** Returns true if the current process has the scene focus */
  virtual bool CanRenderScene() { VRSTUB_NOT_IMPLEMENTED }

  /** Creates a window on the primary monitor to display what is being shown in the headset. */
  virtual void ShowMirrorWindow() { VRSTUB_NOT_IMPLEMENTED }

  /** Closes the mirror window. */
  virtual void HideMirrorWindow() { VRSTUB_NOT_IMPLEMENTED }

  /** Returns true if the mirror window is shown. */
  virtual bool IsMirrorWindowVisible() { VRSTUB_NOT_IMPLEMENTED }

  /** Writes all images that the compositor knows about (including overlays) to a 'screenshots' folder in the SteamVR runtime root. */
  virtual void CompositorDumpImages() { VRSTUB_NOT_IMPLEMENTED }

  /** Let an app know it should be rendering with low resources. */
  virtual bool ShouldAppRenderWithLowResources() { VRSTUB_NOT_IMPLEMENTED }

  /** Override interleaved reprojection logic to force on. */
  virtual void ForceInterleavedReprojectionOn( bool bOverride ) { VRSTUB_NOT_IMPLEMENTED }

  /** Force reconnecting to the compositor process. */
  virtual void ForceReconnectProcess() { VRSTUB_NOT_IMPLEMENTED }

  /** Temporarily suspends rendering (useful for finer control over scene transitions). */
  virtual void SuspendRendering( bool bSuspend ) { VRSTUB_NOT_IMPLEMENTED }

  /** Opens a shared D3D11 texture with the undistorted composited image for each eye.  Use ReleaseMirrorTextureD3D11 when finished
  * instead of calling Release on the resource itself. */
  virtual vr::EVRCompositorError GetMirrorTextureD3D11( vr::EVREye eEye, void *pD3D11DeviceOrResource, void **ppD3D11ShaderResourceView ) { VRSTUB_NOT_IMPLEMENTED }
  virtual void ReleaseMirrorTextureD3D11( void *pD3D11ShaderResourceView ) { VRSTUB_NOT_IMPLEMENTED }

  /** Access to mirror textures from OpenGL. */
  virtual vr::EVRCompositorError GetMirrorTextureGL( vr::EVREye eEye, vr::glUInt_t *pglTextureId, vr::glSharedTextureHandle_t *pglSharedTextureHandle ) { VRSTUB_NOT_IMPLEMENTED }
  virtual bool ReleaseSharedGLTexture( vr::glUInt_t glTextureId, vr::glSharedTextureHandle_t glSharedTextureHandle ) { VRSTUB_NOT_IMPLEMENTED }
  virtual void LockGLSharedTextureForAccess( vr::glSharedTextureHandle_t glSharedTextureHandle ) { VRSTUB_NOT_IMPLEMENTED }
  virtual void UnlockGLSharedTextureForAccess( vr::glSharedTextureHandle_t glSharedTextureHandle ) { VRSTUB_NOT_IMPLEMENTED }

  /** [Vulkan Only]
  * return 0. Otherwise it returns the length of the number of bytes necessary to hold this string including the trailing
  * null.  The string will be a space separated list of-required instance extensions to enable in VkCreateInstance */
  virtual uint32_t GetVulkanInstanceExtensionsRequired( VR_OUT_STRING() char *pchValue, uint32_t unBufferSize ) { VRSTUB_NOT_IMPLEMENTED }

  /** [Vulkan only]
  * return 0. Otherwise it returns the length of the number of bytes necessary to hold this string including the trailing
  * null.  The string will be a space separated list of required device extensions to enable in VkCreateDevice */
  virtual uint32_t GetVulkanDeviceExtensionsRequired( VkPhysicalDevice_T *pPhysicalDevice, VR_OUT_STRING() char *pchValue, uint32_t unBufferSize ) { VRSTUB_NOT_IMPLEMENTED }

  /** [ Vulkan/D3D12 Only ]
  * There are two purposes for SetExplicitTimingMode:
  *	1. To get a more accurate GPU timestamp for when the frame begins in Vulkan/D3D12 applications.
  *	2. (Optional) To avoid having WaitGetPoses access the Vulkan queue so that the queue can be accessed from
  *	another thread while WaitGetPoses is executing.
  *
  * More accurate GPU timestamp for the start of the frame is achieved by the application calling
  * SubmitExplicitTimingData immediately before its first submission to the Vulkan/D3D12 queue.
  * This is more accurate because normally this GPU timestamp is recorded during WaitGetPoses.  In D3D11, 
  * WaitGetPoses queues a GPU timestamp write, but it does not actually get submitted to the GPU until the 
  * application flushes.  By using SubmitExplicitTimingData, the timestamp is recorded at the same place for 
  * Vulkan/D3D12 as it is for D3D11, resulting in a more accurate GPU time measurement for the frame.
  *
  * Avoiding WaitGetPoses accessing the Vulkan queue can be achieved using SetExplicitTimingMode as well.  If this is desired,
  * the application should set the timing mode to Explicit_ApplicationPerformsPostPresentHandoff and *MUST* call PostPresentHandoff
  * itself. If these conditions are met, then WaitGetPoses is guaranteed not to access the queue.  Note that PostPresentHandoff
  * and SubmitExplicitTimingData will access the queue, so only WaitGetPoses becomes safe for accessing the queue from another
  * thread. */
  virtual void SetExplicitTimingMode( EVRCompositorTimingMode eTimingMode ) { VRSTUB_NOT_IMPLEMENTED }

  /** [ Vulkan/D3D12 Only ]
  * Submit explicit timing data.  When SetExplicitTimingMode is true, this must be called immediately before
  * the application's first vkQueueSubmit (Vulkan) or ID3D12CommandQueue::ExecuteCommandLists (D3D12) of each frame.
  * This function will insert a GPU timestamp write just before the application starts its rendering.  This function
  * will perform a vkQueueSubmit on Vulkan so must not be done simultaneously with VkQueue operations on another thread.
  * Returns VRCompositorError_RequestFailed if SetExplicitTimingMode is not enabled. */
  virtual EVRCompositorError SubmitExplicitTimingData() { VRSTUB_NOT_IMPLEMENTED }
};

class VRInputStub: public IVRInput {
public:

  // ---------------  Handle management   --------------- //

  /** Sets the path to the action manifest JSON file that is used by this application. If this information
  * was set on the Steam partner site, calls to this function are ignored. If the Steam partner site
  * setting and the path provided by this call are different, VRInputError_MismatchedActionManifest is returned. 
  * This call must be made before the first call to UpdateActionState or IVRSystem::PollNextEvent. */
  virtual EVRInputError SetActionManifestPath( const char *pchActionManifestPath ) { return VRInputError_None; }

  /** Returns a handle for an action set. This handle is used for all performance-sensitive calls. */
  virtual EVRInputError GetActionSetHandle( const char *pchActionSetName, VRActionSetHandle_t *pHandle ) { *pHandle = (VRActionSetHandle_t)0xFACEFACE; return VRInputError_None; }

  /** Returns a handle for an action. This handle is used for all performance-sensitive calls. */
  virtual EVRInputError GetActionHandle( const char *pchActionName, VRActionHandle_t *pHandle ) { *pHandle = (VRActionHandle_t)0xFACEFACE; return VRInputError_None; }

  /** Returns a handle for any path in the input system. E.g. /user/hand/right */
  virtual EVRInputError GetInputSourceHandle( const char *pchInputSourcePath, VRInputValueHandle_t  *pHandle ) { *pHandle = (VRInputValueHandle_t)0xFACEFACE; return VRInputError_None; }

  // --------------- Reading action state ------------------- //

  /** Reads the current state into all actions. After this call, the results of Get*Action calls 
  * will be the same until the next call to UpdateActionState. */
  virtual EVRInputError UpdateActionState( VR_ARRAY_COUNT( unSetCount ) VRActiveActionSet_t *pSets, uint32_t unSizeOfVRSelectedActionSet_t, uint32_t unSetCount ) { return VRInputError_None; }

  /** Reads the state of a digital action given its handle. This will return VRInputError_WrongType if the type of
  * action is something other than digital */
  virtual EVRInputError GetDigitalActionData( VRActionHandle_t action, InputDigitalActionData_t *pActionData, uint32_t unActionDataSize, VRInputValueHandle_t ulRestrictToDevice );

  /** Reads the state of an analog action given its handle. This will return VRInputError_WrongType if the type of
  * action is something other than analog */
  virtual EVRInputError GetAnalogActionData( VRActionHandle_t action, InputAnalogActionData_t *pActionData, uint32_t unActionDataSize, VRInputValueHandle_t ulRestrictToDevice );

  /** Reads the state of a pose action given its handle. */
  virtual EVRInputError GetPoseActionData( VRActionHandle_t action, ETrackingUniverseOrigin eOrigin, float fPredictedSecondsFromNow, InputPoseActionData_t *pActionData, uint32_t unActionDataSize, VRInputValueHandle_t ulRestrictToDevice );

  /** Reads the state of a skeletal action given its handle. */
  virtual EVRInputError GetSkeletalActionData( VRActionHandle_t action, InputSkeletalActionData_t *pActionData, uint32_t unActionDataSize, VRInputValueHandle_t ulRestrictToDevice ) { VRSTUB_NOT_IMPLEMENTED }

  // --------------- Skeletal Bone Data ------------------- //

  /** Reads the state of the skeletal bone data associated with this action and copies it into the given buffer. */
  virtual EVRInputError GetSkeletalBoneData( VRActionHandle_t action, EVRSkeletalTransformSpace eTransformSpace, EVRSkeletalMotionRange eMotionRange, VR_ARRAY_COUNT( unTransformArrayCount ) VRBoneTransform_t *pTransformArray, uint32_t unTransformArrayCount, VRInputValueHandle_t ulRestrictToDevice ) { VRSTUB_NOT_IMPLEMENTED }

  /** Reads the state of the skeletal bone data in a compressed form that is suitable for
  * sending over the network. The required buffer size will never exceed ( sizeof(VR_BoneTransform_t)*boneCount + 2).
  * Usually the size will be much smaller. */
  virtual EVRInputError GetSkeletalBoneDataCompressed( VRActionHandle_t action, EVRSkeletalTransformSpace eTransformSpace, EVRSkeletalMotionRange eMotionRange, VR_OUT_BUFFER_COUNT( unCompressedSize ) void *pvCompressedData, uint32_t unCompressedSize, uint32_t *punRequiredCompressedSize, VRInputValueHandle_t ulRestrictToDevice ) { VRSTUB_NOT_IMPLEMENTED }

  /** Turns a compressed buffer from GetSkeletalBoneDataCompressed and turns it back into a bone transform array. */
  virtual EVRInputError DecompressSkeletalBoneData( void *pvCompressedBuffer, uint32_t unCompressedBufferSize, EVRSkeletalTransformSpace *peTransformSpace, VR_ARRAY_COUNT( unTransformArrayCount ) VRBoneTransform_t *pTransformArray, uint32_t unTransformArrayCount ) { VRSTUB_NOT_IMPLEMENTED }

  // --------------- Haptics ------------------- //

  /** Triggers a haptic event as described by the specified action */
  virtual EVRInputError TriggerHapticVibrationAction( VRActionHandle_t action, float fStartSecondsFromNow, float fDurationSeconds, float fFrequency, float fAmplitude, VRInputValueHandle_t ulRestrictToDevice ) { VRSTUB_NOT_IMPLEMENTED }

  // --------------- Action Origins ---------------- //

  /** Retrieve origin handles for an action */
  virtual EVRInputError GetActionOrigins( VRActionSetHandle_t actionSetHandle, VRActionHandle_t digitalActionHandle, VR_ARRAY_COUNT( originOutCount ) VRInputValueHandle_t *originsOut, uint32_t originOutCount ) { VRSTUB_NOT_IMPLEMENTED }

  /** Retrieves the name of the origin in the current language */
  virtual EVRInputError GetOriginLocalizedName( VRInputValueHandle_t origin, VR_OUT_STRING() char *pchNameArray, uint32_t unNameArraySize ) { VRSTUB_NOT_IMPLEMENTED }

  /** Retrieves useful information for the origin of this action */
  virtual EVRInputError GetOriginTrackedDeviceInfo( VRInputValueHandle_t origin, InputOriginInfo_t *pOriginInfo, uint32_t unOriginInfoSize ) { VRSTUB_NOT_IMPLEMENTED }

  /** Shows the current binding for the action in-headset */
  virtual EVRInputError ShowActionOrigins( VRActionSetHandle_t actionSetHandle, VRActionHandle_t ulActionHandle ) { VRSTUB_NOT_IMPLEMENTED }

  /** Shows the current binding all the actions in the specified action sets */
  virtual EVRInputError ShowBindingsForActionSet( VR_ARRAY_COUNT( unSetCount ) VRActiveActionSet_t *pSets, uint32_t unSizeOfVRSelectedActionSet_t, uint32_t unSetCount, VRInputValueHandle_t originToHighlight ) { VRSTUB_NOT_IMPLEMENTED }
};


} // stub
} // vr 

#endif /* _H_OpenVRStub */
