#include "Test.h"

#include <glm/gtc/type_ptr.hpp>

#include "Matrix.h"
#include "TTT.h"
#include "Vector.h"
#include "View.h"

static bool IsIdentity(const pymol::TTT& ttt)
{
  auto pre = ttt.getPreTranslation();
  auto rot = ttt.getRotation();
  auto post = ttt.getPostTranslation();
  return pre.x == 0.0f && pre.y == 0.0f && pre.z == 0.0f && post.x == 0.0f &&
         post.y == 0.0f && post.z == 0.0f && rot.x == 0.0f && rot.y == 0.0f &&
         rot.z == 0.0f && rot.w == 1.0f;
}

TEST_CASE("TTT Matrix", "[TTT Matrix]")
{
  pymol::TTT ttt;
  REQUIRE(IsIdentity(ttt));

  pymol::TTT ttt2(glm::vec3(), glm::quat(1.0f, 0.0f, 0.0f, 0.0f), glm::vec3());
  REQUIRE(IsIdentity(ttt2));
}

TEST_CASE("TTT Matrix Reset", "[TTT Matrix]")
{
  pymol::TTT ttt;
  REQUIRE(IsIdentity(ttt));
  ttt.translate(glm::vec3(1.0f, 2.0f, 3.0f));
  REQUIRE(!IsIdentity(ttt));
  ttt.reset();
  REQUIRE(IsIdentity(ttt));
}

TEST_CASE("TTT Matrix Originate & Get Pre Translation", "[TTT Matrix]")
{
  pymol::TTT ttt;
  ttt.translate(glm::vec3(10.0f, 20.0f, 30.0f));
  ttt.originate(glm::vec3(5.0f, 6.0f, 7.0f));
  auto pretranslate = ttt.getPreTranslation();
  REQUIRE(pretranslate.x == -5.0f);
  REQUIRE(pretranslate.y == -6.0f);
  REQUIRE(pretranslate.z == -7.0f);
  auto posttranslate = ttt.getPostTranslation();
  REQUIRE(posttranslate.x == 15.0f);
  REQUIRE(posttranslate.y == 26.0f);
  REQUIRE(posttranslate.z == 37.0f);
  ttt.originate(glm::vec3(100.0f, 100.0f, 100.0f));
  pretranslate = ttt.getPreTranslation();
  REQUIRE(pretranslate.x == -100.0f);
  REQUIRE(pretranslate.y == -100.0f);
  REQUIRE(pretranslate.z == -100.0f);
  posttranslate = ttt.getPostTranslation();
  REQUIRE(posttranslate.x == 110.0f);
  REQUIRE(posttranslate.y == 120.0f);
  REQUIRE(posttranslate.z == 130.0f);
}

TEST_CASE("TTT Matrix Translate & Get Post-Translation", "[TTT Matrix]")
{
  pymol::TTT ttt;
  ttt.translate(glm::vec3(5.0f, 6.0f, 7.0f));
  auto posttranslate = ttt.getPostTranslation();
  REQUIRE(posttranslate.x == 5.0f);
  REQUIRE(posttranslate.y == 6.0f);
  REQUIRE(posttranslate.z == 7.0f);
  float tr[3]{1.0f, 2.0f, 3.0f};
  ttt.translate(tr);
  posttranslate = ttt.getPostTranslation();
  REQUIRE(posttranslate.x == 6.0f);
  REQUIRE(posttranslate.y == 8.0f);
  REQUIRE(posttranslate.z == 10.0f);
}

TEST_CASE("TTT Matrix Get Homogenous Matrix", "[TTT Matrix]")
{
  pymol::TTT ttt;
  auto trans_vec = glm::vec3(5.0f, 6.0f, 7.0f);
  ttt.translate(trans_vec);
  auto hom = ttt.getHomogenousMatrix();
  hom = glm::transpose(hom);
  auto hom_ptr = glm::value_ptr(hom);

  float mat[16];
  identity44f(mat);
  mat[3] += trans_vec[0];
  mat[7] += trans_vec[1];
  mat[11] += trans_vec[2];

  REQUIRE(pymol::almost_equal_n(hom_ptr, 16, mat));
}

TEST_CASE("TTT Matrix Construct from Components", "[TTT Matrix]")
{
  glm::vec3 ori(10.0f, 20.0f, 30.0f);
  auto ori_ptr = glm::value_ptr(ori);
  glm::vec3 tr(5.0f, 15.0f, 25.0f);
  auto ang = glm::radians(90.0f);
  auto axis = glm::vec3(1.0f, 0.0f, 0.0f);
  glm::quat rot = glm::angleAxis(ang, axis);
  auto rot_ptr = glm::value_ptr(rot);
  auto ttt = pymol::TTT(-ori, rot, ori + tr);

  auto n_pre_ttt = -ttt.getPreTranslation();
  auto n_pre_ttt_ptr = glm::value_ptr(n_pre_ttt);
  auto rot_ttt = ttt.getRotation();
  auto rot_ttt_ptr = glm::value_ptr(rot_ttt);
  auto post_ttt = ttt.getPostTranslation();
  auto post_ttt_ptr = glm::value_ptr(post_ttt);

  REQUIRE(pymol::almost_equal_n(ori_ptr, 3, n_pre_ttt_ptr));
  REQUIRE(pymol::almost_equal_n(rot_ptr, 4, rot_ttt_ptr));
  REQUIRE(pymol::almost_equal_n(glm::value_ptr(ori + tr), 3, post_ttt_ptr));
}

TEST_CASE("TTT Matrix Rotate and Get Rotation", "[TTT Matrix]")
{
  pymol::TTT ttt;
  auto axis = glm::vec3(1.0f, 0.0f, 0.0f);
  auto axis2 = glm::vec3(0.0f, 1.0f, 0.0f);
  auto degRot = glm::radians(90.0f);
  ttt.rotate(degRot, axis);
  ttt.rotate(degRot, axis2);
  auto rot = ttt.getRotation();
  REQUIRE(rot[0] == Approx(0.5f));
  REQUIRE(rot[1] == Approx(0.5f));
  REQUIRE(rot[2] == Approx(0.5f));
  REQUIRE(rot[3] == Approx(0.5f));
}

TEST_CASE("TTT Matrix - TTT Matrix Multiply", "[TTT Matrix]")
{
  float tttMat[16];
  identity44f(tttMat);
  tttMat[3] += 3.0f;
  tttMat[7] += 4.0f;
  tttMat[11] += 5.0f;
  float tttMat2[16];
  identity44f(tttMat2);
  tttMat2[3] += 10.0f;
  tttMat2[7] += 20.0f;
  tttMat2[11] += 30.0f;
  float tttMat3[16];
  identity44f(tttMat3);
  combineTTT44f44f(tttMat, tttMat2, tttMat3);

  pymol::TTT ttt;
  ttt.translate(glm::vec3(3.0f, 4.0f, 5.0f));
  pymol::TTT ttt2;
  ttt2.translate(glm::vec3(10.0f, 20.0f, 30.0f));

  auto ttt3 = ttt * ttt2;
  auto posttranslation = ttt3.getPostTranslation();
  REQUIRE(posttranslation[0] == 13.0f);
  REQUIRE(posttranslation[1] == 24.0f);
  REQUIRE(posttranslation[2] == 35.0f);
}

TEST_CASE("TTT Matrix - Transform Position", "[TTT Matrix]")
{
  glm::vec3 ori(100.0f, 100.0f, 100.0f);
  glm::vec3 pos(5.0f, 6.0f, 7.0f);
  float pos2[3]{5.0f, 6.0f, 7.0f};
  float pos2_result[3]{};
  glm::vec3 tr(3.0f, 4.0f, 5.0f);

  pymol::TTT ttt;
  ttt.originate(ori);
  ttt.translate(tr);
  auto newPos = ttt.transform(pos);
  ttt.transform(pos2, pos2_result);


  float tttMat[16];
  identity44f(tttMat);
  float oritttMat[16];
  identity44f(oritttMat);
  oritttMat[3] = ori[0];
  oritttMat[7] = ori[1];
  oritttMat[11] = ori[2];
  oritttMat[12] = -ori[0];
  oritttMat[13] = -ori[1];
  oritttMat[14] = -ori[2];
  float mult[16];
  combineTTT44f44f(tttMat, oritttMat, mult);
  mult[3] = tr[0];
  mult[7] = tr[1];
  mult[11] = tr[2];
  float newPos2[3]{};
  transformTTT44f3f(mult, glm::value_ptr(pos), newPos2);

  REQUIRE(pymol::almost_equal_n(glm::value_ptr(newPos), 3, newPos2));
  REQUIRE(pymol::almost_equal_n(glm::value_ptr(newPos), 3, pos2_result));
}

TEST_CASE("TTT Matrix - Transform Vector", "[TTT Matrix]")
{
  glm::vec3 ori(100.0f, 100.0f, 100.0f);
  glm::vec3 vec(5.0f, 6.0f, 7.0f);
  float vec2[3]{5.0f, 6.0f, 7.0f};
  float vec2_result[3]{};
  glm::vec3 tr(3.0f, 4.0f, 5.0f);

  pymol::TTT ttt;
  ttt.originate(ori);
  ttt.translate(tr);
  auto newVec = ttt.transform_vector(vec);
  ttt.transform_vector(vec2, vec2_result);

  float tttMat[16];
  identity44f(tttMat);
  float oritttMat[16];
  identity44f(oritttMat);
  oritttMat[3] = ori[0];
  oritttMat[7] = ori[1];
  oritttMat[11] = ori[2];
  oritttMat[12] = -ori[0];
  oritttMat[13] = -ori[1];
  oritttMat[14] = -ori[2];
  float mult[16];
  combineTTT44f44f(tttMat, oritttMat, mult);
  mult[3] = tr[0];
  mult[7] = tr[1];
  mult[11] = tr[2];
  float newVec2[3]{};
  transform_normalTTT44f3f(mult, glm::value_ptr(vec), newVec2);

  REQUIRE(pymol::almost_equal_n(glm::value_ptr(newVec), 3, newVec2));
  REQUIRE(pymol::almost_equal_n(glm::value_ptr(newVec), 3, vec2_result));
}

TEST_CASE("TTT Lerp", "[TTT Matrix]")
{
  pymol::TTT ttt1;
  ttt1.translate(glm::vec3(10.0f, 20.0f, 30.0f));
  ttt1.rotate(glm::radians(90.0f), glm::vec3(1.0f, 0.0f, 0.0f));
  ttt1.originate(glm::vec3(100.0f, 200.0f, 300.0f));
  pymol::TTT ttt2;
  ttt2.translate(glm::vec3(20.0f, 30.0f, 40.0f));
  ttt2.rotate(glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
  ttt2.originate(glm::vec3(200.0f, 300.0f, 400.0f));
  auto lerped_ttt = pymol::lerp(ttt1, ttt2, 0.5f);
  auto pre = lerped_ttt.getPreTranslation();
  auto rot = lerped_ttt.getRotation();
  auto post = lerped_ttt.getPostTranslation();

  REQUIRE(pre[0] == Approx(-150.0f));
  REQUIRE(pre[1] == Approx(-250.0f));
  REQUIRE(pre[2] == Approx(-350.0f));
  REQUIRE(pre[0] == Approx(-150.0f));
  REQUIRE(pre[1] == Approx(-250.0f));
  REQUIRE(pre[2] == Approx(-350.0f));
  REQUIRE(rot[0] == Approx(0.40825f));
  REQUIRE(rot[1] == Approx(0.40825f));
  REQUIRE(rot[2] == Approx(0.0f));
  REQUIRE(rot[3] == Approx(0.40825f * 2.0f));
  REQUIRE(post[0] == Approx( 165.0f));
  REQUIRE(post[1] == Approx( 275.0f));
  REQUIRE(post[2] == Approx( 385.0f));
}

TEST_CASE("TTT Serialize To", "[TTT Matrix]")
{
  glm::vec3 ori(100.0f, 100.0f, 100.0f);
  glm::vec3 vec(5.0f, 6.0f, 7.0f);
  glm::vec3 tr(3.0f, 4.0f, 5.0f);

  pymol::TTT ttt;
  ttt.originate(ori);
  ttt.translate(tr);

  auto ttt_result = pymol::TTT::as_pymol_2_legacy(ttt);
  auto ttt_ptr = glm::value_ptr(ttt_result);

  float tttMat[16];
  identity44f(tttMat);
  tttMat[3] = ori[0];
  tttMat[7] = ori[1];
  tttMat[11] = ori[2];
  tttMat[12] = -ori[0];
  tttMat[13] = -ori[1];
  tttMat[14] = -ori[2];
  tttMat[3] += tr[0];
  tttMat[7] += tr[1];
  tttMat[11] += tr[2]; // translate

  REQUIRE(pymol::almost_equal_n(ttt_ptr, 16, tttMat));
}

TEST_CASE("TTT Serialize From", "[TTT Matrix]")
{
  glm::vec3 ori(100.0f, 100.0f, 100.0f);
  glm::vec3 vec(5.0f, 6.0f, 7.0f);
  glm::vec3 tr(3.0f, 4.0f, 5.0f);

  float ttt[16];
  identity44f(ttt);
  ttt[3] = ori[0];
  ttt[7] = ori[1];
  ttt[11] = ori[2];
  ttt[12] = -ori[0];
  ttt[13] = -ori[1];
  ttt[14] = -ori[2];
  ttt[3] += tr[0];
  ttt[7] += tr[1];
  ttt[11] += tr[2]; // translate

  auto newTTT = pymol::TTT::from_pymol_2_legacy(ttt);
  auto result = pymol::TTT::as_pymol_2_legacy(newTTT);
  auto ttt_ptr = glm::value_ptr(result);

  REQUIRE(pymol::almost_equal_n(ttt_ptr, 16, ttt));
}

TEST_CASE("TTT Get Rotation About", "[TTT Matrix]")
{
  glm::vec3 ori(5.0f, 10.0f, 15.f);
  glm::vec3 dir(1.0f, 0.0f, 0.0f);
  float angle = glm::radians(90.0f);
  auto ttt = pymol::TTT::rotation_about_with_origin(angle, dir, ori);
  auto ttt_leg = pymol::TTT::as_pymol_2_legacy(ttt);
  auto ttt_ptr = glm::value_ptr(ttt_leg);

  float old_ttt[16];
  get_rotation_about3f3fTTTf(
      angle, glm::value_ptr(dir), glm::value_ptr(ori), old_ttt);

  REQUIRE(pymol::almost_equal_n(ttt_ptr, 16, old_ttt));
}

TEST_CASE("TTT to SceneView", "[TTT Matrix]")
{
  glm::vec3 ori(100.0f, 100.0f, 100.0f);
  glm::vec3 pos(5.0f, 6.0f, 7.0f);
  glm::vec3 tr(3.0f, 4.0f, 5.0f);
  float ang = glm::radians(45.0f);
  auto axis = glm::vec3(1.0f, 0.0f, 0.0f);

  pymol::TTT ttt;
  ttt.originate(ori);
  ttt.translate(tr);
  ttt.rotate(ang, axis);

  float tttMat[16];
  identity44f(tttMat);
  tttMat[3] = ori[0];
  tttMat[7] = ori[1];
  tttMat[11] = ori[2];
  tttMat[12] = -ori[0];
  tttMat[13] = -ori[1];
  tttMat[14] = -ori[2];
  tttMat[3] += tr[0];
  tttMat[7] += tr[1];
  tttMat[11] += tr[2]; // translate
  glm::mat4 rotMat = glm::transpose(glm::mat4_cast(glm::angleAxis(ang, axis)));
  rotMat[0][3] = ori[0];
  rotMat[1][3] = ori[1];
  rotMat[2][3] = ori[2];
  rotMat[3][0] = -ori[0];
  rotMat[3][1] = -ori[1];
  rotMat[3][2] = -ori[2];
  float res_ttt[16];
  combineTTT44f44f(tttMat, glm::value_ptr(rotMat), res_ttt);

  auto elem1 = pymol::TTT::to_view_elem(ttt);

  CViewElem elem2{};
  TTTToViewElem(res_ttt, &elem2);

  REQUIRE(pymol::almost_equal_n(elem1.matrix, 16, elem2.matrix));
  REQUIRE(pymol::almost_equal_n(elem1.pre, 3, elem2.pre));
  REQUIRE(pymol::almost_equal_n(elem1.post, 3, elem2.post));
}

TEST_CASE("SceneView to TTT", "[TTT Matrix]")
{
  glm::vec3 ori(100.0f, 100.0f, 100.0f);
  glm::vec3 pos(5.0f, 6.0f, 7.0f);
  glm::vec3 tr(3.0f, 4.0f, 5.0f);
  float ang = glm::radians(45.0f);
  auto axis = glm::vec3(1.0f, 0.0f, 0.0f);

  pymol::TTT ttt;
  ttt.originate(ori);
  ttt.translate(tr);
  ttt.rotate(ang, axis);

  auto elem_fwd = pymol::TTT::to_view_elem(ttt);
  auto elem_rev = pymol::TTT::from_view_elem(elem_fwd);
  auto elem_fwd2 = pymol::TTT::to_view_elem(elem_rev);
  REQUIRE(pymol::almost_equal_n(elem_fwd.matrix, 16, elem_fwd2.matrix));
  REQUIRE(pymol::almost_equal_n(elem_fwd.pre, 3, elem_fwd2.pre));
  REQUIRE(pymol::almost_equal_n(elem_fwd.post, 3, elem_fwd2.post));
}
