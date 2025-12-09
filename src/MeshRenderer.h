#pragma once

#include <vector>

#if defined(_WIN32) || defined(_WIN64)
#include <Windows.h>
#endif
#include <GL/gl.h>

// 頂点座標を格納する構造体です。
struct Vertex {
    GLfloat x;
    GLfloat y;
    GLfloat z;
};

// 面を頂点インデックスで格納する構造体です。
struct Face {
    GLuint v1;
    GLuint v2;
    GLuint v3;
};

// 線分を頂点インデックスで格納する構造体です。
struct Line {
    GLuint v1;
    GLuint v2;
};

// 頂点に設定する色データを格納する構造体です。
struct Color {
    GLfloat r;
    GLfloat g;
    GLfloat b;
    GLfloat a;
};

// 頂点に対応する法線ベクトルを格納する構造体です。
struct Normal {
    GLfloat x;
    GLfloat y;
    GLfloat z;
};

/**
 * メッシュおよび線分を描画するためのクラスです。
 * */
class MeshRenderer
{
private:
    /* 定数定義 */

    // 矢印メッシュの円錐部分の全体に占める長さの割合です。
    const double ARROW_CONE_LENGTH_RATIO = 0.1;
    // 矢印メッシュの軸部分の全体に占める長さの割合です。
    const double ARROW_SHAFT_LENGTH_RATIO = 0.9;

    // 矢印メッシュの円錐部分の半径です。
    const GLfloat ARROW_CONE_RADIUS = 0.03;
    // 矢印メッシュの軸部分の半径です。
    const GLfloat ARROW_SHAFT_RADIUS = 0.01;
    
    // 矢印の直線部分と円錐部分のメッシュの分割数です。
    // NOTE: 描画速度を優先し、最小限にしています。
    const int ARROW_SLICES = 3;

    // 円柱の半径を補正するための係数です。
    const double CYLINDER_RADIUS_CORRECT_RATIO = 1e-2;


private:
    /* メッシュ用の配列 */

    // メッシュの頂点データを格納する配列です。
    std::vector<Vertex> mesh_vertex_list;
    // メッシュの面データを格納する配列です。
    std::vector<Face> mesh_face_list;
    // 頂点毎の色データを格納する配列です。
    std::vector<Color> mesh_color_list;
    // 頂点毎の法線ベクトルを格納する配列です。
    std::vector<Normal> mesh_normal_list;

    // mesh_vertex_list 配列に格納されている有効な頂点の数です。
    int mesh_vertex_count = 0;
    // mesh_face_list 配列に格納されている有効な面の数です。
    int mesh_face_count = 0;


    /* 線分用の配列 */

    // 線分の頂点データを格納する配列です。
    std::vector<Vertex> line_vertex_list;
    // 線分データを格納する配列です。
    std::vector<Line> line_list;

    // line_vertex_list 配列に格納されている有効な頂点の数です。
    int line_vertex_count = 0;
    // line_list 配列に格納されている有効な線分の数です。
    int line_count = 0;


public:

    /**
     * メッシュ、線分の登録を全て解除します。
     * */
    void clear();


    /**
     * 球体のメッシュを描画対象として追加します。
     * @param x         球体の中心 X 座標です。
     * @param y         球体の中心 Y 座標です。
     * @param z         球体の中心 Z 座標です。
     * @param r         球体の半径です。
     * @param color     球体の色 (RGBA 4ch 分のデータを持つ配列) です。
     * @param division  球体のメッシュの分割数。大きいほどなめらかになります。
     * */
    void appendSphere(
        GLfloat x,
        GLfloat y,
        GLfloat z,
        GLfloat r,
        const GLfloat* color,
        int division);

    /**
     * 矢印のメッシュを描画対象として追加します。
     * @param x        矢印の開始 X 座標です。
     * @param y        矢印の開始 Y 座標です。
     * @param z        矢印の開始 Z 座標です。
     * @param vx       矢印の向きを表すベクトルの X です。
     * @param vy       矢印の向きを表すベクトルの Y です。
     * @param vz       矢印の向きを表すベクトルの Z です。
     * @param scale    表示倍率。矢印の太さ以外に影響します。
     * */
    void appendArrowVector(
        GLfloat x, GLfloat y, GLfloat z,
        GLfloat vx, GLfloat vy, GLfloat vz,
        float scale = 1.0);

    /**
     * 矢印のメッシュを描画対象として追加します。
     * @param x0       矢印の開始 X 座標です。
     * @param y0       矢印の開始 Y 座標です。
     * @param z0       矢印の開始 Z 座標です。
     * @param x1       矢印の終端 X 座標です。
     * @param y1       矢印の終端 Y 座標です。
     * @param z1       矢印の終端 Z 座標です。
     * */
    void appendArrowPosition(
        GLfloat x0, GLfloat y0, GLfloat z0,
        GLfloat x1, GLfloat y1, GLfloat z1);

    /**
     * 三角形を描画対象として追加します。
     * @param x0       頂点 0 の X 座標です。
     * @param y0       頂点 0 の Y 座標です。
     * @param z0       頂点 0 の Z 座標です。
     * @param x1       頂点 1 の X 座標です。
     * @param y1       頂点 1 の Y 座標です。
     * @param z1       頂点 1 の Z 座標です。
     * @param x2       頂点 2 の X 座標です。
     * @param y2       頂点 2 の Y 座標です。
     * @param z2       頂点 2 の Z 座標です。
     * */
    void appendTriangle(
        GLfloat x0, GLfloat y0, GLfloat z0,
        GLfloat x1, GLfloat y1, GLfloat z1,
        GLfloat x2, GLfloat y2, GLfloat z2);

    /**
     * 線分を 1 つ配列に追加します。
     * @param x0       頂点 0 の X 座標です。
     * @param y0       頂点 0 の Y 座標です。
     * @param z0       頂点 0 の Z 座標です。
     * @param x1       頂点 1 の X 座標です。
     * @param y1       頂点 1 の Y 座標です。
     * @param z1       頂点 1 の Z 座標です。
     * */
    void appendLine(
        GLfloat x0, GLfloat y0, GLfloat z0,
        GLfloat x1, GLfloat y1, GLfloat z1);

    /**
     * 円柱のメッシュを配列に追加します。
     * @param x0           円柱の開始 X 座標です。
     * @param y0           円柱の開始 Y 座標です。
     * @param z0           円柱の開始 Z 座標です。
     * @param x1           円柱の終了 X 座標です。
     * @param y1           円柱の終了 Y 座標です。
     * @param z1           円柱の終了 Z 座標です。
     * @param radius       円柱の半径です。(この値に補正係数を適用した値が半径となります)
     * @param div_number   円柱の分割数です。
     * */
    void appendCylinder(
        GLfloat x0,
        GLfloat y0,
        GLfloat z0,
        GLfloat x1,
        GLfloat y1,
        GLfloat z1,
        float radius,
        int div_number = 8);


    /**
     * 配列に登録されているメッシュを描画します。
     * */
    void drawMeshesWithColor();

    /**
     * 配列に登録されているメッシュを描画します。
     * */
    void drawMeshesWithNormal();

    /**
     * 配列に登録されているメッシュを描画します。
     * */
    void drawMeshes();

    /**
     * 配列に登録されているメッシュをワイヤーフレームとして描画します。
     * */
    void drawWireframe();

    /**
     * 配列に登録されている線分を描画します。
     * @param line_width       描画する線の太さです。
     * */
    void drawLines(GLfloat line_width = 1.0);

};

