#include "myglutdisplay.h"


// 別のファイル (util.cpp) で宣言されている関数です。クォータニオンから対応する回転行列を計算します。
void qua2rot(float r[], float q[]);


/* 2024/03 CJS 追加ここから */
#include "MeshRenderer.h"

/* 定数定義 */

// 球メッシュの最小分割数です。3 以上を指定します。
const int MIN_DIV_NUMBER = 3;
// 原点からの距離がこの値付近の原子から、メッシュが荒くなります。
const double DISTANCE_THRESHOLD = 11.0;

/**
 * クォータニオンとスクロール量から画面中心に来る座標を計算します。
 * @param [in] quaternion    クォータニオンです。
 * @param [in] position      スクロール量です。
 * @param [out] cx           中心座標 X を受け取る変数です。
 * @param [out] cy           中心座標 Y を受け取る変数です。
 * @param [out] cz           中心座標 Z を受け取る変数です。
 * */
void calculateDisplayCenter(
    float quaternion[],
    float position[],
    float& cx,
    float& cy,
    float& cz)
{
    // 逆クォータニオンおよび逆回転行列を納める配列です。
    float inverse_quaternion[4];
    float inverse_rotation_matrix[16];

    // 逆クォータニオンを計算します。
    float sum = quaternion[0] * quaternion[0] + quaternion[1] * quaternion[1] + quaternion[2] * quaternion[2] + quaternion[3] * quaternion[3];
    if (!sum) sum = 1;
    inverse_quaternion[0] = quaternion[0] / sum;
    inverse_quaternion[1] = -quaternion[1] / sum;
    inverse_quaternion[2] = -quaternion[2] / sum;
    inverse_quaternion[3] = -quaternion[3] / sum;

    // 逆クォータニオンから逆回転行列を作成します。
    qua2rot(inverse_rotation_matrix, inverse_quaternion);

    // スクロール量を反転させた値を逆回転の適用対象座標とします。
    float bx = -position[0];
    float by = -position[1];
    float bz = 1;

    // 逆回転行列を適用し、画面の中央に来るべき座標を決定します。
    cx = inverse_rotation_matrix[0] * bx + inverse_rotation_matrix[4] * by + inverse_rotation_matrix[8] * bz;
    cy = inverse_rotation_matrix[1] * bx + inverse_rotation_matrix[5] * by + inverse_rotation_matrix[9] * bz;
    cz = inverse_rotation_matrix[2] * bx + inverse_rotation_matrix[6] * by + inverse_rotation_matrix[10] * bz;
}

/**
 * 原子の座標と表示されている画像の中心からの距離に基づき、最終的な球メッシュの粒度を決定します。
 * @param [in] divide_number    基準となるメッシュの粒度です。
 * @param [in] cx               中心座標の X です。
 * @param [in] cy               中心座標の Y です。
 * @param [in] cz               中心座標の Z です。
 * @param [in] x                原子座標の X です。
 * @param [in] y                原子座標の Y です。
 * @param [in] z                原子座標の Z です。
 * @param [in] display_scale    表示倍率です。
 * @return                      球メッシュの粒度を返します。
 * */
inline int calcSphereDivide(
    int divide_number,
    float cx, float cy, float cz,
    float x, float y, float z,
    float display_scale)
{
    float distance = std::sqrt((cx - x) * (cx - x) + (cy - y) * (cy - y) + (cz - z) * (cz - z));

    // シグモイド関数の性質を指定する係数です。
    const double sigmoid_k = 10.0;
    int atom_div_number = divide_number / (1.0 + exp(sigmoid_k * (distance * display_scale - DISTANCE_THRESHOLD)));

    return (atom_div_number >= MIN_DIV_NUMBER) ? atom_div_number : MIN_DIV_NUMBER;
}

/* 2024/03 CJS 追加ここまで */

/* 2024/11 CJS 追加ここから */

/**
 * ヒット機能がオンの時は原子の描画を行い、
 * ヒット機能がオフの時は描画クラスへの原子の登録が行います。
 * @param enable_atom_click     ヒット機能のオン/オフを表すフラグです。0: ヒット機能オフ 1: ヒット機能オン
 * @param xx                    原子の X 座標です。
 * @param yy                    原子の Y 座標です。
 * @param zz                    原子の Z 座標です。
 * @param rad                   原子の半径です。
 * @param atom_color            原子の描画色です。
 * @param scale                 表示倍率です。
 * @param atom_div_count        原子のメッシュ分割数です。値が大きいほどメッシュの頂点数が増え、球に近くなります。
 * @param pname                 メッシュ表面のマテリアル パラメーターです。
 * @param renderer              ヒット機能がオフの時の原子の登録先となる描画クラスです。
 * */
void drawAtom(
    int enable_atom_click,
    float xx, float yy, float zz, float rad,
    const GLfloat* atom_color, float scale, int atom_div_count,
    GLenum pname, MeshRenderer& renderer)
{
    if (enable_atom_click) {
        glMaterialfv(GL_FRONT_AND_BACK, pname, atom_color);
        glTranslated(xx * scale, yy * scale, zz * scale);
        glutSolidSphere(rad * scale, atom_div_count, atom_div_count);
        glTranslated(-xx * scale, -yy * scale, -zz * scale);
    }
    else {
        renderer.appendSphere(xx * scale, yy * scale, zz * scale, rad * scale, atom_color, atom_div_count);
    }
}

/* 2024/03 CJS 追加ここまで */

void drawCell(MeshRenderer& renderer);
void drawYZPlane(MeshRenderer& renderer, float x);
void glDrawAxisd(MeshRenderer& renderer, float length);
void cnt_switch_visibility_by_wall();


void myGlutDisplay(void)
{
    // 描画クラス。頂点データ等を収めるメモリの再確保を防ぐため、static としています。
    static MeshRenderer renderer;

    // 球の分割レベルを表示倍率に合わせて変更します。
    // 倍率が小さいときはポリゴン数を減らし、
    // 倍率が高いときは segments で設定されている値に近づけます。
    GLuint sphere_div_number = (GLuint)(segments * (scl / 3));
    if (sphere_div_number < MIN_DIV_NUMBER) sphere_div_number = MIN_DIV_NUMBER;
    if (sphere_div_number > segments) sphere_div_number = segments;

    // 画面の中心に来る座標を算出します。
    float cx;
    float cy;
    float cz;
    calculateDisplayCenter(cq, obj_pos, cx, cy, cz);

    /* 2024/03 CJS 追加ここまで */


    float xx1, yy1, zz1, xx2, yy2, zz2, xx3, yy3, zz3, xx4, yy4, zz4;
    glClearColor(.9f, .9f, .9f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    //  glMatrixMode( GL_PROJECTION );
    //  glMatrixMode( GL_MODELVIEW );
    glPushMatrix();
    glTranslatef(obj_pos[0], obj_pos[1], -obj_pos[2]);
    glMultMatrixf(rotmat);

    float rad = (float)radius / 100.0;

    // draw cell
    if (show_cell) {
        drawCell(renderer);
    }
    // draw axis
    if (show_axis) {
        glDisable(GL_LIGHTING);
        glDrawAxisd(renderer, 0.5);
    }

    if (ievec) {
        if (ievec_num < 1) { ievec_num = 1; }
        if (ievec_num > MAXMODE) { ievec_num = MAXMODE; }
        //glDisable( GL_LIGHTING );  glColor3d(1.0, 0.0, 1.0);
        glEnable(GL_LIGHTING);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, green);

        // 描画クラスに登録されているデータをクリアします。2024/03 CJS
        renderer.clear();

        for (int i = 1; i <= atom.natom; i++) {
            float xx1 = atom.rx[i] * 1e9; float yy1 = atom.ry[i] * 1e9; float zz1 = atom.rz[i] * 1e9;
            float xx2 = atom.evecx[i][ievec_num - 1] * evec_len;
            float yy2 = atom.evecy[i][ievec_num - 1] * evec_len;
            float zz2 = atom.evecz[i][ievec_num - 1] * evec_len;

            // 矢印メッシュを描画クラスに登録します。2024/03 CJS
            renderer.appendArrowVector(xx1, yy1, zz1, xx2, yy2, zz2, scl);
        }

        // 描画クラスに登録したメッシュを描画します。2024/03 CJS
        renderer.drawMeshes();

        //glEnable( GL_LIGHTING );
    }

    // draw elements
    if (atom.QC) {
        glDisable(GL_LIGHTING);
        glColor3d(0.0, 0.5, 0.5);

        // 描画クラスに登録されているデータをクリアします。2024/03 CJS
        renderer.clear();

        for (int i = 1; i <= atom.nelem; i++) {
            mk_helmat(i, atom.rx, atom.ry, atom.rz, mat);
            int ix0 = 0; int iy0 = 0; int iz0 = 0;
            if (atom.elem_v_rep[i][1] >= 4) { ix0 = 1; }
            else { ix0 = 0; }
            if ((atom.elem_v_rep[i][1] % 4 == 2) || (atom.elem_v_rep[i][1] % 4 == 3)) { iy0 = 1; }
            else { iy0 = 0; }
            if (atom.elem_v_rep[i][1] % 2 == 1) { iz0 = 1; }
            else { iz0 = 0; }
            xx1 = atom.rx[atom.elem_v[i][1]] * 1e9; yy1 = atom.ry[atom.elem_v[i][1]] * 1e9; zz1 = atom.rz[atom.elem_v[i][1]] * 1e9;
            xx1 = xx1 + (cell.hmat[0][0] * ix0 + cell.hmat[0][1] * iy0 + cell.hmat[0][2] * iz0) * 1e9;
            yy1 = yy1 + (cell.hmat[1][0] * ix0 + cell.hmat[1][1] * iy0 + cell.hmat[1][2] * iz0) * 1e9;
            zz1 = zz1 + (cell.hmat[2][0] * ix0 + cell.hmat[2][1] * iy0 + cell.hmat[2][2] * iz0) * 1e9;
            xx2 = xx1 + mat[0][0] * 1e9; yy2 = yy1 + mat[1][0] * 1e9; zz2 = zz1 + mat[2][0] * 1e9;
            xx3 = xx1 + mat[0][1] * 1e9; yy3 = yy1 + mat[1][1] * 1e9; zz3 = zz1 + mat[2][1] * 1e9;
            xx4 = xx1 + mat[0][2] * 1e9; yy4 = yy1 + mat[1][2] * 1e9; zz4 = zz1 + mat[2][2] * 1e9;

            // 線分を描画クラスに登録します。2024/03 CJS
            renderer.appendLine(xx1 * scl, yy1 * scl, zz1 * scl, xx2 * scl, yy2 * scl, zz2 * scl);
            renderer.appendLine(xx2 * scl, yy2 * scl, zz2 * scl, xx3 * scl, yy3 * scl, zz3 * scl);
            renderer.appendLine(xx1 * scl, yy1 * scl, zz1 * scl, xx3 * scl, yy3 * scl, zz3 * scl);

            renderer.appendLine(xx1 * scl, yy1 * scl, zz1 * scl, xx4 * scl, yy4 * scl, zz4 * scl);
            renderer.appendLine(xx2 * scl, yy2 * scl, zz2 * scl, xx4 * scl, yy4 * scl, zz4 * scl);
            renderer.appendLine(xx3 * scl, yy3 * scl, zz3 * scl, xx4 * scl, yy4 * scl, zz4 * scl);
        }

        // 描画クラスに登録した線分を描画します。2024/03 CJS
        renderer.drawLines();

        glEnable(GL_LIGHTING);
    }

    // draw walls (for CNT)
    if (show_cnt_wall) {
        glEnable(GL_LIGHTING);
        glColor3d(0.0, 0.5, 0.5);
        float cube3 = (float)cell.hmat[2][2] * 1e9;

        // 描画クラスに登録されているデータをクリアします。2024/03 CJS
        renderer.clear();

        for (int i = 1; i <= atom.nwall; i++) {
            int i0 = atom.wall_v[i][0]; int i1 = atom.wall_v[i][1]; int i2 = atom.wall_v[i][2];
            if ((i0 != 0) && (i1 != 0) && (i2 != 0)) {
                xx1 = atom.rx[abs(i0)] * 1e9; yy1 = atom.ry[abs(i0)] * 1e9; zz1 = atom.rz[abs(i0)] * 1e9;
                xx2 = atom.rx[abs(i1)] * 1e9; yy2 = atom.ry[abs(i1)] * 1e9; zz2 = atom.rz[abs(i1)] * 1e9;
                xx3 = atom.rx[abs(i2)] * 1e9; yy3 = atom.ry[abs(i2)] * 1e9; zz3 = atom.rz[abs(i2)] * 1e9;
                if (zz2 > zz1 + cube3 / 2) { zz2 = zz2 - cube3; } if (zz2 < zz1 - cube3 / 2) { zz2 = zz2 + cube3; }
                if (zz3 > zz1 + cube3 / 2) { zz3 = zz3 - cube3; } if (zz3 < zz1 - cube3 / 2) { zz3 = zz3 + cube3; }
                if (cube3 < 3.2) {
                    if (i1 < 0) { if (zz2 < cube3 / 2) { zz2 += cube3; } else { zz2 -= cube3; } }
                    if (i2 < 0) { if (zz3 < cube3 / 2) { zz3 += cube3; } else { zz3 -= cube3; } }
                }

                // ポリゴンを 1 つ描画クラスに登録します。2024/03 CJS
                renderer.appendTriangle(
                    xx1* scl, yy1* scl, zz1* scl,
                    xx2* scl, yy2* scl, zz2* scl,
                    xx3* scl, yy3* scl, zz3* scl);
            }
        }

        // 描画クラスに登録したメッシュを描画します。2024/03 CJS
        renderer.drawMeshes();


        glDisable(GL_LIGHTING);
        glColor3d(0.5, 0.5, 0.5);

        // 描画クラスに登録されているデータをクリアします。2024/03 CJS
        renderer.clear();

        for (int i = 1; i <= atom.nwall; i++) {
            int i0 = atom.wall_v[i][0]; int i1 = atom.wall_v[i][1]; int i2 = atom.wall_v[i][2];
            if ((i0 != 0) && (i1 != 0) && (i2 != 0)) {
                xx1 = atom.rx[abs(i0)] * 1e9; yy1 = atom.ry[abs(i0)] * 1e9; zz1 = atom.rz[abs(i0)] * 1e9;
                xx2 = atom.rx[abs(i1)] * 1e9; yy2 = atom.ry[abs(i1)] * 1e9; zz2 = atom.rz[abs(i1)] * 1e9;
                xx3 = atom.rx[abs(i2)] * 1e9; yy3 = atom.ry[abs(i2)] * 1e9; zz3 = atom.rz[abs(i2)] * 1e9;
                if (zz2 > zz1 + cube3 / 2) { zz2 = zz2 - cube3; } if (zz2 < zz1 - cube3 / 2) { zz2 = zz2 + cube3; }
                if (zz3 > zz1 + cube3 / 2) { zz3 = zz3 - cube3; } if (zz3 < zz1 - cube3 / 2) { zz3 = zz3 + cube3; }
                if (cube3 < 3.2) {
                    if (i1 < 0) { if (zz2 < cube3 / 2) { zz2 += cube3; } else { zz2 -= cube3; } }
                    if (i2 < 0) { if (zz3 < cube3 / 2) { zz3 += cube3; } else { zz3 -= cube3; } }
                }

                // 線分を描画クラスに登録します。2024/03 CJS
                renderer.appendLine(xx1* scl, yy1* scl, zz1* scl, xx2* scl, yy2* scl, zz2* scl);
                renderer.appendLine(xx2* scl, yy2* scl, zz2* scl, xx3* scl, yy3* scl, zz3* scl);
                renderer.appendLine(xx3* scl, yy3* scl, zz3* scl, xx1* scl, yy1* scl, zz1* scl);

                if (show_cnt_wallv) {
                    float xx0, yy0, zz0, vx, vy, vz, vscale = vscl;
                    xx0 = (xx1 + xx2 + xx3) / 3; yy0 = (yy1 + yy2 + yy3) / 3; zz0 = (zz1 + zz2 + zz3) / 3;
                    vx = atom.wall_nvec[i][0]; vy = atom.wall_nvec[i][1]; vz = atom.wall_nvec[i][2];
                    vx = xx0 + vx * vscale; vy = yy0 + vy * vscale; vz = zz0 + vz * vscale;

                    // 矢印のメッシュを描画クラスに登録します。2024/03 CJS
                    renderer.appendArrowPosition(xx0* scl, yy0* scl, zz0* scl, vx* scl, vy* scl, vz* scl);
                }
            }
        }

        // 描画クラスに登録したメッシュと線分を描画します。2024/03 CJS
        renderer.drawMeshes();
        renderer.drawLines();

    }

    if ((show_cnt_wall_num != 0) && (show_cnt_wall_num <= atom.nwall)) {
        glEnable(GL_LIGHTING);
        glColor3d(1.0, 1.0, 1.0);
        float cube3 = (float)cell.hmat[2][2] * 1e9;
        int i = show_cnt_wall_num;
        int i0 = atom.wall_v[i][0]; int i1 = atom.wall_v[i][1]; int i2 = atom.wall_v[i][2];
        if ((i0 != 0) && (i1 != 0) && (i2 != 0)) {
            // 描画クラスに登録されているデータをクリアします。2024/03 CJS
            renderer.clear();

            xx1 = atom.rx[abs(i0)] * 1e9; yy1 = atom.ry[abs(i0)] * 1e9; zz1 = atom.rz[abs(i0)] * 1e9;
            xx2 = atom.rx[abs(i1)] * 1e9; yy2 = atom.ry[abs(i1)] * 1e9; zz2 = atom.rz[abs(i1)] * 1e9;
            xx3 = atom.rx[abs(i2)] * 1e9; yy3 = atom.ry[abs(i2)] * 1e9; zz3 = atom.rz[abs(i2)] * 1e9;
            if (zz2 > zz1 + cube3 / 2) { zz2 = zz2 - cube3; } if (zz2 < zz1 - cube3 / 2) { zz2 = zz2 + cube3; }
            if (zz3 > zz1 + cube3 / 2) { zz3 = zz3 - cube3; } if (zz3 < zz1 - cube3 / 2) { zz3 = zz3 + cube3; }
            if (cube3 < 3.2) {
                if (i1 < 0) { if (zz2 < cube3 / 2) { zz2 += cube3; } else { zz2 -= cube3; } }
                if (i2 < 0) { if (zz3 < cube3 / 2) { zz3 += cube3; } else { zz3 -= cube3; } }
            }

            // ポリゴンを 1 つ描画クラスに登録します。2024/03 CJS
            renderer.appendTriangle(
                xx1 * scl, yy1 * scl, zz1 * scl,
                xx2 * scl, yy2 * scl, zz2 * scl,
                xx3 * scl, yy3 * scl, zz3 * scl);

            float xx0, yy0, zz0, vx, vy, vz, vscale = vscl;
            xx0 = (xx1 + xx2 + xx3) / 3; yy0 = (yy1 + yy2 + yy3) / 3; zz0 = (zz1 + zz2 + zz3) / 3;
            vx = atom.wall_nvec[i][0]; vy = atom.wall_nvec[i][1]; vz = atom.wall_nvec[i][2];
            vx = xx0 + vx * vscale; vy = yy0 + vy * vscale; vz = zz0 + vz * vscale;

            // 矢印のメッシュを描画クラスに登録します。2024/03 CJS
            renderer.appendArrowPosition(xx0* scl, yy0* scl, zz0* scl, vx* scl, vy* scl, vz* scl);

            // 登録済みのメッシュを描画します。2024/03 CJS
            renderer.drawMeshes();
        }

        glDisable(GL_LIGHTING);
        glColor3d(0.5, 0.5, 0.5);
        //int i = show_cnt_wall_num;
        //int i0 = atom.wall_v[i][0]; int i1 = atom.wall_v[i][1]; int i2 = atom.wall_v[i][2];
        if ((i0 != 0) && (i1 != 0) && (i2 != 0)) {
            // 描画クラスに登録されているデータをクリアします。2024/03 CJS
            renderer.clear();

            xx1 = atom.rx[abs(i0)] * 1e9; yy1 = atom.ry[abs(i0)] * 1e9; zz1 = atom.rz[abs(i0)] * 1e9;
            xx2 = atom.rx[abs(i1)] * 1e9; yy2 = atom.ry[abs(i1)] * 1e9; zz2 = atom.rz[abs(i1)] * 1e9;
            xx3 = atom.rx[abs(i2)] * 1e9; yy3 = atom.ry[abs(i2)] * 1e9; zz3 = atom.rz[abs(i2)] * 1e9;
            if (zz2 > zz1 + cube3 / 2) { zz2 = zz2 - cube3; } if (zz2 < zz1 - cube3 / 2) { zz2 = zz2 + cube3; }
            if (zz3 > zz1 + cube3 / 2) { zz3 = zz3 - cube3; } if (zz3 < zz1 - cube3 / 2) { zz3 = zz3 + cube3; }
            if (cube3 < 3.2) {
                if (i1 < 0) { if (zz2 < cube3 / 2) { zz2 += cube3; } else { zz2 -= cube3; } }
                if (i2 < 0) { if (zz3 < cube3 / 2) { zz3 += cube3; } else { zz3 -= cube3; } }
            }

            // 線分を描画クラスに登録します。2024/03 CJS
            renderer.appendLine(xx1* scl, yy1* scl, zz1* scl, xx2* scl, yy2* scl, zz2* scl);
            renderer.appendLine(xx2* scl, yy2* scl, zz2* scl, xx3* scl, yy3* scl, zz3* scl);
            renderer.appendLine(xx3* scl, yy3* scl, zz3* scl, xx1* scl, yy1* scl, zz1* scl);

            // 描画クラスに登録した線分を描画します。2024/03 CJS
            renderer.drawLines();
        }
    }
    // draw auxiliary constraint
    if (draw_aux) {
        if (yzplane_punch) {
            drawYZPlane(renderer, yzplane_punch_d / 2);
            drawYZPlane(renderer , -yzplane_punch_d / 2);
        }
    }
    // draw ring (for CNT)
    if (show_cnt_ring) {
        GLUquadricObj* ring = gluNewQuadric();
        gluQuadricDrawStyle(ring, GLU_LINE);
        float cx = (float)cell.hmat[0][0];
        float cy = (float)cell.hmat[1][1];
        float cz = (float)cell.hmat[2][2];
        float rad = (float)cnt_ring_radius * ang;
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, red);
        glPushMatrix();
        glTranslated(cx / 2 * 1e9 * scl, cy / 2 * 1e9 * scl, 0);
        //glTranslated( 0,0,0);
        gluCylinder(ring, rad * 1e9 * scl, rad * 1e9 * scl, cz * 1e9 * scl, 50, 20);
        gluDeleteQuadric(ring);
        glPopMatrix();
    }

    if (draw_bond) {
        if (bond_display_type == 0) {
            if (bond_display_color == 0) {
                glColor3f(red[0], red[1], red[2]);
            }
            else if (bond_display_color == 1) {
                glColor3f(blue[0], blue[1], blue[2]);
            }
            else if (bond_display_color == 2) {
                glColor3f(green[0], green[1], green[2]);
            }
            else if (bond_display_color == 3) {
                glColor3f(yellow[0], yellow[1], yellow[2]);
            }
            else if (bond_display_color == 4) {
                glColor3f(purple[0], purple[1], purple[2]);
            }
            else if (bond_display_color == 5) {
                glColor3f(gray[0], gray[1], gray[2]);
            }
            else if (bond_display_color == 6) {
                glColor3f(black[0], black[1], black[2]);
            }
            else if (bond_display_color == 7) {
                glColor3f(white[0], white[1], white[2]);
            }
        }
        else {
            if (bond_display_color == 0) {
                glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, red);
            }
            else if (bond_display_color == 1) {
                glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, blue);
            }
            else if (bond_display_color == 2) {
                glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, green);
            }
            else if (bond_display_color == 3) {
                glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, yellow);
            }
            else if (bond_display_color == 4) {
                glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, purple);
            }
            else if (bond_display_color == 5) {
                glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, gray);
            }
            else if (bond_display_color == 6) {
                glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, black);
            }
            else if (bond_display_color == 7) {
                glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, white);
            }
        }
        glDisable(GL_LIGHTING);
        //glColor3d(0.5, 0.5, 0.5);
        float cube1 = (float)cell.hmat[0][0] * 1e9;
        float cube2 = (float)cell.hmat[1][1] * 1e9;
        float cube3 = (float)cell.hmat[2][2] * 1e9;
        float cube1x = (float)cell.hmat[0][0] * 1e9;
        float cube1y = (float)cell.hmat[1][0] * 1e9;
        float cube1z = (float)cell.hmat[2][0] * 1e9;
        float cube2x = (float)cell.hmat[0][1] * 1e9;
        float cube2y = (float)cell.hmat[1][1] * 1e9;
        float cube2z = (float)cell.hmat[2][1] * 1e9;
        float cube3x = (float)cell.hmat[0][2] * 1e9;
        float cube3y = (float)cell.hmat[1][2] * 1e9;
        float cube3z = (float)cell.hmat[2][2] * 1e9;
        if (book.algo == 3) {
            //for (int i=1; i<=bre.np; i++) {

            // 描画クラスに登録されているデータをクリアします。2024/03 CJS
            renderer.clear();

            for (int i = 1; i <= atom.natom; i++) {
                if (!atom.visible[i]) { continue; }
                if (atom.nneighbor[i] > MAXNEIGHBOR) { continue; }
                for (int j = 0; j < atom.nneighbor[i]; j++) {
                    int jn = atom.neighbor[i][j];
                    if (jn > atom.natom) { continue; }
                    xx1 = atom.rx[i] * 1e9; yy1 = atom.ry[i] * 1e9; zz1 = atom.rz[i] * 1e9;
                    xx2 = atom.rx[jn] * 1e9; yy2 = atom.ry[jn] * 1e9; zz2 = atom.rz[jn] * 1e9;
                    int drawflg = 1;
                    if (cell.pbcx) {
                        if (xx2 > xx1 + cube1 / 1.5) {
                            xx2 = xx2 - cube1; if (draw_bond_pbc == 0) { drawflg = 0; }
                        }
                        if (xx2 < xx1 - cube1 / 1.5) {
                            xx2 = xx2 + cube1; if (draw_bond_pbc == 0) { drawflg = 0; }
                        }
                    }
                    if (cell.pbcy) {
                        if (yy2 > yy1 + cube2 / 1.5) {
                            yy2 = yy2 - cube2; if (draw_bond_pbc == 0) { drawflg = 0; }
                        }
                        if (yy2 < yy1 - cube2 / 1.5) {
                            yy2 = yy2 + cube2; if (draw_bond_pbc == 0) { drawflg = 0; }
                        }
                    }
                    if (cell.pbcz) {
                        if (zz2 > zz1 + cube3 / 1.5) {
                            zz2 = zz2 - cube3; if (draw_bond_pbc == 0) { drawflg = 0; }
                        }
                        if (zz2 < zz1 - cube3 / 1.5) {
                            zz2 = zz2 + cube3; if (draw_bond_pbc == 0) { drawflg = 0; }
                        }
                    }

                    if (drawflg == 1) {
                        if (bond_display_type == 0) { // Simple lines
                            // 線分を描画クラスに登録します。2024/03 CJS
                            renderer.appendLine(xx1* scl, yy1* scl, zz1* scl, xx2* scl, yy2* scl, zz2* scl);

                        }
                        else { // Cylinder with lighting
                            // 円柱を描画クラスに登録します。2024/03 CJS
                            renderer.appendCylinder(
                                xx1* scl, yy1* scl, zz1* scl, xx2* scl, yy2* scl, zz2* scl,
                                bondthickness* scl / 5.);
                        }
                    }
                }
            }

            // 描画クラスに登録した線分と円柱を描画します。円柱はライティングを有効にします。2024/03 CJS 追加ここから
            if (bond_display_type == 0) {
                renderer.drawLines(bondthickness);
            }
            else {
                glEnable(GL_LIGHTING);
                renderer.drawMeshesWithNormal();
                glDisable(GL_LIGHTING);
            }
            // 2024/03 CJS 追加ここまで
        }
        else if (book.algo == 1) { // Bug fix YU 20140508
            //std::cout<<"[myglutdisplay.cpp]"<<std::endl;
            //std::cout<<"Setting bond using BK algo #1"<<std::endl;
            //std::cout<<"Not supported yet."<<std::endl;
            float blang = bondlength * 1.0e-10;

            // 描画クラスに登録されているデータをクリアします。2024/03 CJS
            renderer.clear();

            for (int i = 1; i <= atom.natom; i++) {
                if (atom.nneighbor[i] > MAXNEIGHBOR) { continue; }
                for (int j = 0; j < atom.nneighbor[i]; j++) {
                    int jn = atom.neighbor[i][j];
                    xx1 = atom.rx[i] * 1e9; yy1 = atom.ry[i] * 1e9; zz1 = atom.rz[i] * 1e9;
                    xx2 = atom.rx[jn] * 1e9; yy2 = atom.ry[jn] * 1e9; zz2 = atom.rz[jn] * 1e9;
                    xx3 = atom.rx[jn] - atom.rx[i]; yy3 = atom.ry[jn] - atom.ry[i]; zz3 = atom.rz[jn] - atom.rz[i];
                    xx4 = cell.hinmat[0][0] * xx3 + cell.hinmat[0][1] * yy3 + cell.hinmat[0][2] * zz3;
                    yy4 = cell.hinmat[1][0] * xx3 + cell.hinmat[1][1] * yy3 + cell.hinmat[1][2] * zz3;
                    zz4 = cell.hinmat[2][0] * xx3 + cell.hinmat[2][1] * yy3 + cell.hinmat[2][2] * zz3;
                    int drawflg = 1; int chk = 0;
                    if (cell.pbcx) {
                        if (xx4 > 0.5) { xx2 -= cube1x; yy2 -= cube1y; zz2 -= cube1z; chk = 1; }
                        if (xx4 < -0.5) { xx2 += cube1x; yy2 += cube1y; zz2 += cube1z; chk = 1; }
                        //if (xx3 > blang) {xx2-=cube1x;yy2-=cube1y;zz2-=cube1z;chk=1;}
                        //if (xx3 <-blang) {xx2+=cube1x;yy2+=cube1y;zz2+=cube1z;chk=1;}
                    }
                    else {
                        //if (xx4> 0.5) {drawflg=0;} if (xx4<-0.5) {drawflg=0;}
                        if (fabs(xx3) > blang) { drawflg = 0; }
                    }
                    if (cell.pbcy) {
                        if (yy4 > 0.5) { xx2 -= cube2x; yy2 -= cube2y; zz2 -= cube2z; chk = 1; }
                        if (yy4 < -0.5) { xx2 += cube2x; yy2 += cube2y; zz2 += cube2z; chk = 1; }
                        //if (yy3 > blang) {xx2-=cube2x;yy2-=cube2y;zz2-=cube2z;chk=1;}
                        //if (yy3 <-blang) {xx2+=cube2x;yy2+=cube2y;zz2+=cube2z;chk=1;}
                    }
                    else {
                        //if (yy4> 0.5) {drawflg=0;} if (yy4<-0.5) {drawflg=0;}
                        if (fabs(yy3) > blang) { drawflg = 0; }
                    }
                    if (cell.pbcz) {
                        if (zz4 > 0.5) { xx2 -= cube3x; yy2 -= cube3y; zz2 -= cube3z; chk = 1; }
                        if (zz4 < -0.5) { xx2 += cube3x; yy2 += cube3y; zz2 += cube3z; chk = 1; }
                        //if (zz3 > blang) {xx2-=cube3x;yy2-=cube3y;zz2-=cube3z;chk=1;}
                        //if (zz3 <-blang) {xx2+=cube3x;yy2+=cube3y;zz2+=cube3z;chk=1;}
                    }
                    else {
                        //if (zz4> 0.5) {drawflg=0;} if (zz4<-0.5) {drawflg=0;}
                        if (fabs(zz3) > blang) { drawflg = 0; }
                    }
                    if ((chk == 1) && (draw_bond_pbc == 0)) { drawflg = 0; }

                    if (drawflg == 1) {
                        if (bond_display_type == 0) { // Simple lines
                            // 線分を描画クラスに登録します。2024/03 CJS
                            renderer.appendLine(xx1* scl, yy1* scl, zz1* scl, xx2* scl, yy2* scl, zz2* scl);
                        }
                        else { // Cylinder with lighting
                            // 円柱を描画クラスに登録します。2024/03 CJS
                            renderer.appendCylinder(xx1* scl, yy1* scl, zz1* scl, xx2* scl, yy2* scl, zz2* scl, bondthickness* scl / 5.);
                        }
                        //
                    }
                }
            }

            // 描画クラスに登録した線分と円柱を描画します。円柱はライティングを有効にします。2024/03 CJS 追加ここから
            if (bond_display_type == 0) {
                renderer.drawLines(bondthickness);
            }
            else {
                glEnable(GL_LIGHTING);
                renderer.drawMeshesWithNormal();
                glDisable(GL_LIGHTING);
            }
            // 2024/03 CJS 追加ここまで

        }
    } // end of draw_bond

    glEnable(GL_LIGHTING);

    // draw force
    if (draw_force) {
        // 描画クラスに登録されているデータをクリアします。2024/03 CJS
        renderer.clear();

        glEnable(GL_LIGHTING);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, green2);
        for (int i = 1; i <= atom.natom; i++) {
            float xx0, yy0, zz0, vx, vy, vz, vscale;
            xx0 = atom.rx[i] * 1e9; yy0 = atom.ry[i] * 1e9; zz0 = atom.rz[i] * 1e9;
            vx = atom.fx[i];     vy = atom.fy[i];     vz = atom.fz[i];
            double vlen = sqrt(vx * vx + vy * vy + vz * vz);
            vscale = 0.5 / (0.1 * eV / ang) * vscl_force; //0.1eV/A -> length 0.5
            //if (vlen>0.1*eV/ang) { vscale=0.5/vlen*vscl_force; } //Force over 0.1eV/A -> length 0.5
            //vscale=0.5/atom.Fmax(); //Fmax -> length 0.5
            //vx = xx0 + vx*vscale; vy = yy0 + vy*vscale; vz = zz0 + vz*vscale; 
            vx = vx * vscale; vy = vy * vscale; vz = vz * vscale;
            //glDrawArrowd(xx0*scl,yy0*scl,zz0*scl,vx*scl,vy*scl,vz*scl);
            
            // 矢印メッシュを描画クラスに登録します。2024/03 CJS
            renderer.appendArrowVector(xx0, yy0, zz0, vx, vy, vz, scl);
        }

        // 描画クラスに登録したメッシュを描画します。2024/03 CJS
        renderer.drawMeshes();
    } // end of draw force

    if (draw_load) {
        // 描画クラスに登録されているデータをクリアします。2024/03 CJS
        renderer.clear();

        glEnable(GL_LIGHTING);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, blue2);
        for (int i = 1; i <= atom.natom; i++) {
            float xx0, yy0, zz0, vx, vy, vz, vscale;
            xx0 = atom.rx[i] * 1e9; yy0 = atom.ry[i] * 1e9; zz0 = atom.rz[i] * 1e9;
            vx = atom.fx_l[i];   vy = atom.fy_l[i];   vz = atom.fz_l[i];
            double vlen = sqrt(vx * vx + vy * vy + vz * vz);
            vscale = 0.5 / (0.1 * eV / ang) * vscl_force; //0.1eV/A -> length 0.5
            //if (vlen>0.1*eV/ang) { vscale=0.5/vlen*vscl_force; } //Force over 0.1eV/A -> length 0.5
            //vscale=0.5/atom.Fmax(); //Fmax -> length 0.5
            //vx = xx0 + vx*vscale; vy = yy0 + vy*vscale; vz = zz0 + vz*vscale; 
            vx = vx * vscale; vy = vy * vscale; vz = vz * vscale;
            //glDrawArrowd(xx0*scl,yy0*scl,zz0*scl,vx*scl,vy*scl,vz*scl);
            
            // 矢印メッシュを描画クラスに登録します。2024/03 CJS 
            renderer.appendArrowVector(xx0, yy0, zz0, vx, vy, vz, scl);
        }

        // 描画クラスに登録したメッシュを描画します。2024/03 CJS
        renderer.drawMeshes();
    } // end of draw force
    glPopMatrix();

    // draw atoms
    //  float rad = (float)radius/100.0;
    if (mdmotion) {
        change_atom_color();
    }
    if ((createconfig_mode > 0)||(editatom_mode > 0)) {
        if ((iatom_pick > 0) && (iatom_pick <= atom.natom)) {
            for (int i = 1; i <= atom.natom; i++) {
                if (i == iatom_pick) {
                    memcpy(color[i], green, sizeof(GLfloat) * 4);
                }
                else {
                    memcpy(color[i], color0[i], sizeof(GLfloat) * 4);
                }
            }
        }
    }
    //if (mode_cnt_corrugation&&(!show_cnt_all)) {
    if (mode_cnt_corrugation) { //YU2025 (this 'if' bifurcation may need checking for CNT mode)
    for (int i = 1; i <= atom.natom; i++) { atom.visible[i] = true; }
    if (!show_cnt_all) {
        cnt_switch_visibility_by_wall();
    }
    }

    // 描画クラスに登録されているデータをクリアします。2024/03 CJS
    renderer.clear();

    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);
    int counter = 0;

    for (int i = 1; i <= atom.natom; i++) {
        glNewList(objects + i, GL_COMPILE);

        // 2024/08 CJS 原子描画前に移動行列を適用します。
        // NOTE: あたり判定で利用する glNewList の処理に回転を含める必要があるため、
        //       原子ごとに回転、移動行列を適用し直します。
        glPushMatrix();
        glTranslatef(obj_pos[0], obj_pos[1], -obj_pos[2]);
        glMultMatrixf(rotmat);

        float xx = atom.rx[i] * 1e9;
        float yy = atom.ry[i] * 1e9;
        float zz = atom.rz[i] * 1e9;

        /* 2024/03 CJS 以下、分岐構造は残したまま、描画処理を appendSphere 関数に差し替えました。*/

        // 原子の画面中心からの距離に基づき、分割数を最終決定します。
        int atom_div_count = calcSphereDivide(sphere_div_number, cx, cy, cz, xx, yy, zz, scl);

        // 原子を表す球メッシュを描画します。
        // 原子のヒット機能がオンの場合は直接描画を行い、
        // 原子のヒット機能がオフの場合は、描画クラスに原子のメッシュを登録した後、一括して描画を行います。2024/11 CJS
        if (atom.QC) {
            int iflg = 0;
            if (edit_elem_mode > 0) {
                for (int ii = 0; ii < 4; ii++) { if ((select_atom[ii] == i) && (select_atom_repidx[ii] == 0)) { iflg = 1; } }
            }
            if (show_only_elem) {
                if (atom.repatom[i] == 0) {
                    // 原子を描画します。 2024/11 CJS
                    drawAtom(enable_atom_click, xx, yy, zz, rad, white, scl, atom_div_count, GL_DIFFUSE, renderer);
                }
            }
            else if (atom.repatom[i] == 1) {
                if (iflg == 0) {
                    // 原子を描画します。 2024/11 CJS
                    drawAtom(enable_atom_click, xx, yy, zz, rad, red, scl, atom_div_count, GL_DIFFUSE, renderer);
                }
                else {
                    // 原子を描画します。 2024/11 CJS
                    drawAtom(enable_atom_click, xx, yy, zz, rad, color[i], scl, atom_div_count, GL_DIFFUSE, renderer);
                }
            }
            else {
                if (iflg == 0) {
                    // 原子を描画します。 2024/11 CJS
                    drawAtom(enable_atom_click, xx, yy, zz, rad, white, scl, atom_div_count, GL_DIFFUSE, renderer);
                }
                else {
                    // 原子を描画します。 2024/11 CJS
                    drawAtom(enable_atom_click, xx, yy, zz, rad, color[i], scl, atom_div_count, GL_DIFFUSE, renderer);
                }
            }
        }
        else { //if not atom.QC
            if (atom.visible[i]) {
                // 原子を描画します。 2024/11 CJS
                drawAtom(enable_atom_click, xx, yy, zz, rad, color[i], scl, atom_div_count, GL_AMBIENT_AND_DIFFUSE, renderer);
            }
        }

        glPopMatrix();
        glEndList();
        glCallList(objects + i);

        /* 2024/03 CJS 更新ここまで。*/
    } // end of loop for real atoms

    // 2024/03 CJS 追加ここから。
    // 2024/08 CJS 移動行列の適用箇所を移動しました。
    glPushMatrix();
    glTranslatef(obj_pos[0], obj_pos[1], -obj_pos[2]);
    glMultMatrixf(rotmat);

    // 原子のヒット機能がオフの場合、描画対象として登録された原子のメッシュを一括で描画します。2024/11 CJS
    if (!enable_atom_click) renderer.drawMeshesWithColor();

    // 変換行列を開放します。
    glPopMatrix();
    // 2024/03 CJS 追加ここまで。

    ibase = objects + atom.natom;
    if (draw_replica > 0) { // replica atoms

        // 描画クラスに登録されているデータをクリアします。2024/03 CJS
        renderer.clear();

        icnt = 0; int ix0, iy0, iz0;
        for (int rbit = 1; rbit <= 7; rbit++) {
            if (rbit >= 4) { ix0 = 1; }
            else { ix0 = 0; }
            if ((rbit % 4 == 2) || (rbit % 4 == 3)) { iy0 = 1; }
            else { iy0 = 0; }
            if (rbit % 2 == 1) { iz0 = 1; }
            else { iz0 = 0; }
            if ((draw_replica == 2) && (iy0 == 1)) { continue; }
            for (int i = 1; i <= atom.natom; i++) {
                float xx = (atom.rx[i] + cell.hmat[0][0] * ix0) * 1e9;
                float yy = (atom.ry[i] + cell.hmat[1][1] * iy0) * 1e9;
                float zz = (atom.rz[i] + cell.hmat[2][2] * iz0) * 1e9;
                if ((xx / 1e9 / cell.hmat[0][0] < 1.1) && (yy / 1e9 / cell.hmat[1][1] < 1.1) && (zz / 1e9 / cell.hmat[2][2] < 1.1) && (icnt < atom.natom * 2)) {
                    icnt++;
                    iatom[icnt] = i; repidx[icnt] = rbit;

                    // 原子のヒット機能がオンかオフかで描画方法を変更します。 2024/11 CJS
                    if (enable_atom_click) {
                        glNewList(ibase + icnt, GL_COMPILE); //printf("glNewList (2)\n");
                        glPushMatrix();
                        glTranslatef(obj_pos[0], obj_pos[1], -obj_pos[2]);
                        glMultMatrixf(rotmat);

                        int atom_div_count = calcSphereDivide(sphere_div_number, cx, cy, cz, xx, yy, zz, scl);
                        glutWireSphere(rad* scl, atom_div_count, atom_div_count);

                        // 変換行列を開放します。
                        glPopMatrix();
                        glEndList();
                        glCallList(ibase + icnt);
                    }
                    else {
                        renderer.appendSphere(
                            xx * scl, yy * scl, zz * scl, rad * scl, color[atom.natom + icnt],
                            calcSphereDivide(sphere_div_number, cx, cy, cz, xx, yy, zz, scl));
                    }
                }
            }
        }

        // 2024/03 CJS 追加ここから。
        // 登録されたスフィアの描画を実行します。
        // 変換行列を登録します (オリジナル コードの移植です)。
        glPushMatrix();
        glTranslatef(obj_pos[0], obj_pos[1], -obj_pos[2]);
        glMultMatrixf(rotmat);

        // 原子のヒット機能がオフの場合、描画対象として登録された原子のメッシュを一括で描画します。2024/11 CJS
        if (!enable_atom_click) renderer.drawWireframe();

        // 変換行列を開放します。
        glPopMatrix();
        // 2024/03 CJS 追加ここまで。


        // Draw lines between selected atoms
        glPushMatrix();
        glTranslatef(obj_pos[0], obj_pos[1], -obj_pos[2]);
        glMultMatrixf(rotmat);
        glDisable(GL_LIGHTING); glColor3d(0.0, 0.0, 0.0);

        // 描画クラスに登録されているデータをクリアします。2024/03 CJS
        renderer.clear();

        for (int ii = 0; ii < 3; ii++) {
            for (int jj = ii + 1; jj < 4; jj++) {
                int si = select_atom[ii]; int sj = select_atom[jj];
                if ((si != 0) && (sj != 0)) {
                    int ix0, iy0, iz0, rbit;
                    rbit = select_atom_repidx[ii];
                    if (rbit >= 4) { ix0 = 1; }
                    else { ix0 = 0; }
                    if ((rbit % 4 == 2) || (rbit % 4 == 3)) { iy0 = 1; }
                    else { iy0 = 0; }
                    if (rbit % 2 == 1) { iz0 = 1; }
                    else { iz0 = 0; }
                    float xx1 = (atom.rx[si] + cell.hmat[0][0] * ix0) * 1e9;
                    float yy1 = (atom.ry[si] + cell.hmat[1][1] * iy0) * 1e9;
                    float zz1 = (atom.rz[si] + cell.hmat[2][2] * iz0) * 1e9;
                    rbit = select_atom_repidx[jj];
                    if (rbit >= 4) { ix0 = 1; }
                    else { ix0 = 0; }
                    if ((rbit % 4 == 2) || (rbit % 4 == 3)) { iy0 = 1; }
                    else { iy0 = 0; }
                    if (rbit % 2 == 1) { iz0 = 1; }
                    else { iz0 = 0; }
                    float xx2 = (atom.rx[sj] + cell.hmat[0][0] * ix0) * 1e9;
                    float yy2 = (atom.ry[sj] + cell.hmat[1][1] * iy0) * 1e9;
                    float zz2 = (atom.rz[sj] + cell.hmat[2][2] * iz0) * 1e9;

                    // 線分を描画クラスに登録します。2024/03 CJS 変更
                    renderer.appendLine(xx1* scl, yy1* scl, zz1* scl, xx2* scl, yy2* scl, zz2* scl);
                }
            }
        }

        // 描画クラスに登録された線分を描画します。2024/03 CJS
        renderer.drawLines();

        glEnable(GL_LIGHTING);
        glPopMatrix();
    } // end of draw_replica
    else if (replica_range > 0.0) { // if not draw_replica but replica_range > 0

        // 描画クラスに登録されているデータをクリアします。2024/03 CJS
        renderer.clear();

        int ix0, iy0, iz0;
        for (int i = 1; i <= atom.natom; i++) {
            float qx = cell.hinmat[0][0] * atom.rx[i] + cell.hinmat[0][1] * atom.ry[i] + cell.hinmat[0][2] * atom.rz[i];
            float qy = cell.hinmat[1][0] * atom.rx[i] + cell.hinmat[1][1] * atom.ry[i] + cell.hinmat[1][2] * atom.rz[i];
            float qz = cell.hinmat[2][0] * atom.rx[i] + cell.hinmat[2][1] * atom.ry[i] + cell.hinmat[2][2] * atom.rz[i];
            ix0 = 0; iy0 = 0; iz0 = 0;
            if (qx < replica_range) ix0 = 1;
            if (qx > 1. - replica_range) ix0 = -1;
            if (qy < replica_range) iy0 = 1;
            if (qy > 1. - replica_range) iy0 = -1;
            if (qz < replica_range) iz0 = 1;
            if (qz > 1. - replica_range) iz0 = -1;
            if (!cell.pbcx) ix0 = 0;
            if (!cell.pbcy) iy0 = 0;
            if (!cell.pbcz) iz0 = 0;
            if ((ix0 == 0) && (iy0 == 0) && (iz0 == 0)) { continue; }
            for (int ii = 0; ii <= abs(ix0); ii++) {
                for (int jj = 0; jj <= abs(iy0); jj++) {
                    for (int kk = 0; kk <= abs(iz0); kk++) {
                        float xx = atom.rx[i]
                            + cell.hmat[0][0] * (float)ix0 * abs(ii)
                            + cell.hmat[0][1] * (float)iy0 * abs(jj)
                            + cell.hmat[0][2] * (float)iz0 * abs(kk);
                        float yy = atom.ry[i]
                            + cell.hmat[1][0] * (float)ix0 * abs(ii)
                            + cell.hmat[1][1] * (float)iy0 * abs(jj)
                            + cell.hmat[1][2] * (float)iz0 * abs(kk);
                        float zz = atom.rz[i]
                            + cell.hmat[2][0] * (float)ix0 * abs(ii)
                            + cell.hmat[2][1] * (float)iy0 * abs(jj)
                            + cell.hmat[2][2] * (float)iz0 * abs(kk);
                        xx *= 1.0e9; yy *= 1.0e9; zz *= 1.0e9;

                        // 原子のメッシュを描画クラスに登録します。2024/03 CJS
                        renderer.appendSphere(
                            xx* scl, yy* scl, zz* scl, rad* scl, color[i],
                            calcSphereDivide(sphere_div_number, cx, cy, cz, xx, yy, zz, scl));
                    }
                }
            }
        }

        // 2024/03 CJS 追加ここから
        // 登録されたスフィアの描画を実行します。
        // 変換行列を登録します (オリジナル コードの移植です)。
        glPushMatrix();
        glTranslatef(obj_pos[0], obj_pos[1], -obj_pos[2]);
        glMultMatrixf(rotmat);

        // 描画を実行します。
        //renderer.drawWireframe();
        if (enable_atom_click) renderer.drawWireframe(); //YU2024.12.20

        // 変換行列を開放します。
        glPopMatrix();
        // 2024/03 CJS 追加ここまで。

    } // end of replica_range > 0.0

    glLoadIdentity();
    gluLookAt(0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
}

void drawCell(MeshRenderer& renderer)
{
    // この関数は描画クラスを用いた実装に変更しました。2024/03 CJS

    // 描画クラスの登録済みメッシュをクリアします。
    renderer.clear();

    float xx1, yy1, zz1, xx2, yy2, zz2;

    glDisable(GL_LIGHTING);
    glColor3d(0.0, 0.0, 0.0);

    // 線分の座標を計算し、描画クラスに登録します。

    xx1 = 0.0; yy1 = 0.0; zz1 = 0.0;
    xx2 = cell.hmat[0][0] * 1e9; yy2 = cell.hmat[1][0] * 1e9; zz2 = cell.hmat[2][0] * 1e9;
    renderer.appendLine(xx1 * scl, yy1 * scl, zz1 * scl, xx2 * scl, yy2 * scl, zz2 * scl);
    xx1 = 0.0; yy1 = 0.0; zz1 = 0.0;
    xx2 = cell.hmat[0][1] * 1e9; yy2 = cell.hmat[1][1] * 1e9; zz2 = cell.hmat[2][1] * 1e9;
    renderer.appendLine(xx1 * scl, yy1 * scl, zz1 * scl, xx2 * scl, yy2 * scl, zz2 * scl);
    xx1 = 0.0; yy1 = 0.0; zz1 = 0.0;
    xx2 = cell.hmat[0][2] * 1e9; yy2 = cell.hmat[1][2] * 1e9; zz2 = cell.hmat[2][2] * 1e9;
    renderer.appendLine(xx1 * scl, yy1 * scl, zz1 * scl, xx2 * scl, yy2 * scl, zz2 * scl);

    xx1 = cell.hmat[0][0] * 1e9; yy1 = cell.hmat[1][0] * 1e9; zz1 = cell.hmat[2][0] * 1e9;
    xx2 = xx1 + cell.hmat[0][1] * 1e9; yy2 = yy1 + cell.hmat[1][1] * 1e9; zz2 = zz1 + cell.hmat[2][1] * 1e9;
    renderer.appendLine(xx1 * scl, yy1 * scl, zz1 * scl, xx2 * scl, yy2 * scl, zz2 * scl);
    xx2 = xx1 + cell.hmat[0][2] * 1e9; yy2 = yy1 + cell.hmat[1][2] * 1e9; zz2 = zz1 + cell.hmat[2][2] * 1e9;
    renderer.appendLine(xx1 * scl, yy1 * scl, zz1 * scl, xx2 * scl, yy2 * scl, zz2 * scl);
    xx1 = xx2; yy1 = yy2; zz1 = zz2;
    xx2 = xx1 + cell.hmat[0][1] * 1e9; yy2 = yy1 + cell.hmat[1][1] * 1e9; zz2 = zz1 + cell.hmat[2][1] * 1e9;
    renderer.appendLine(xx1 * scl, yy1 * scl, zz1 * scl, xx2 * scl, yy2 * scl, zz2 * scl);
    xx1 = xx2; yy1 = yy2; zz1 = zz2;
    xx2 = xx1 - cell.hmat[0][2] * 1e9; yy2 = yy1 - cell.hmat[1][2] * 1e9; zz2 = zz1 - cell.hmat[2][2] * 1e9;
    renderer.appendLine(xx1 * scl, yy1 * scl, zz1 * scl, xx2 * scl, yy2 * scl, zz2 * scl);

    xx1 = cell.hmat[0][1] * 1e9; yy1 = cell.hmat[1][1] * 1e9; zz1 = cell.hmat[2][1] * 1e9;
    xx2 = xx1 + cell.hmat[0][0] * 1e9; yy2 = yy1 + cell.hmat[1][0] * 1e9; zz2 = zz1 + cell.hmat[2][0] * 1e9;
    renderer.appendLine(xx1 * scl, yy1 * scl, zz1 * scl, xx2 * scl, yy2 * scl, zz2 * scl);
    xx2 = xx1 + cell.hmat[0][2] * 1e9; yy2 = yy1 + cell.hmat[1][2] * 1e9; zz2 = zz1 + cell.hmat[2][2] * 1e9;
    renderer.appendLine(xx1 * scl, yy1 * scl, zz1 * scl, xx2 * scl, yy2 * scl, zz2 * scl);
    xx1 = xx2; yy1 = yy2; zz1 = zz2;
    xx2 = xx1 + cell.hmat[0][0] * 1e9; yy2 = yy1 + cell.hmat[1][0] * 1e9; zz2 = zz1 + cell.hmat[2][0] * 1e9;
    renderer.appendLine(xx1 * scl, yy1 * scl, zz1 * scl, xx2 * scl, yy2 * scl, zz2 * scl);

    xx1 = cell.hmat[0][2] * 1e9; yy1 = cell.hmat[1][2] * 1e9; zz1 = cell.hmat[2][2] * 1e9;
    xx2 = xx1 + cell.hmat[0][0] * 1e9; yy2 = yy1 + cell.hmat[1][0] * 1e9; zz2 = zz1 + cell.hmat[2][0] * 1e9;
    renderer.appendLine(xx1 * scl, yy1 * scl, zz1 * scl, xx2 * scl, yy2 * scl, zz2 * scl);
    xx2 = xx1 + cell.hmat[0][1] * 1e9; yy2 = yy1 + cell.hmat[1][1] * 1e9; zz2 = zz1 + cell.hmat[2][1] * 1e9;
    renderer.appendLine(xx1 * scl, yy1 * scl, zz1 * scl, xx2 * scl, yy2 * scl, zz2 * scl);

    // 描画クラスに登録された線分を描画します。
    renderer.drawLines();
    glEnable(GL_LIGHTING);
}

void drawYZPlane(MeshRenderer& renderer, float x)
{
    // この関数は描画クラスを用いた実装に変更しました。2024/03 CJS

    // 描画クラスの登録済みメッシュをクリアします。
    renderer.clear();

    float xx1, yy1, zz1, xx2, yy2, zz2;
    glDisable(GL_LIGHTING);
    glColor3d(1.0, 0.0, 1.0);

    // 線分の座標を計算し、描画クラスに登録します。

    for (int i = 0; i <= 10; i++) {
        xx1 = cell.hmat[0][0] / 2 + x * ang; yy1 = 0.0; zz1 = cell.hmat[2][2] / 10 * (float)i;
        xx2 = cell.hmat[0][0] / 2 + x * ang; yy2 = cell.hmat[1][1]; zz2 = cell.hmat[2][2] / 10 * (float)i;
        xx1 *= 1e9; yy1 *= 1e9; zz1 *= 1e9; xx2 *= 1e9; yy2 *= 1e9; zz2 *= 1e9;
        renderer.appendLine(xx1 * scl, yy1 * scl, zz1 * scl, xx2 * scl, yy2 * scl, zz2 * scl);
    }
    for (int i = 0; i <= 10; i++) {
        xx1 = cell.hmat[0][0] / 2 + x * ang; yy1 = cell.hmat[1][1] / 10 * (float)i; zz1 = 0.0;
        xx2 = cell.hmat[0][0] / 2 + x * ang; yy2 = cell.hmat[1][1] / 10 * (float)i; zz2 = cell.hmat[2][2];
        xx1 *= 1e9; yy1 *= 1e9; zz1 *= 1e9; xx2 *= 1e9; yy2 *= 1e9; zz2 *= 1e9;
        renderer.appendLine(xx1 * scl, yy1 * scl, zz1 * scl, xx2 * scl, yy2 * scl, zz2 * scl);
    }

    // 描画クラスに登録された線分を描画します。
    renderer.drawLines();
    glEnable(GL_LIGHTING);
}

void glDrawAxisd(MeshRenderer& renderer, float length) // Draw Coordinate Axis
{
    GLUquadricObj* arrows[3];

    // Draw X-axis
    glColor3ub(255, 0, 0);
    // 線分を描画します。2024/03 CJS
    renderer.clear();
    renderer.appendLine(-length * scl * 0.5, -0.1 * scl, -0.1 * scl, length * scl, -0.1 * scl, -0.1 * scl);
    renderer.drawLines();

    glPushMatrix();
    arrows[0] = gluNewQuadric();
    gluQuadricDrawStyle(arrows[0], GLU_FILL);
    glTranslated(length * scl, -0.1 * scl, -0.1 * scl);
    glRotated(90.0f, 0, 1, 0);
    gluCylinder(arrows[0], length / 30 * scl, 0.0f, length / 5 * scl, 8, 8);
    gluDeleteQuadric(arrows[0]);
    glPopMatrix();

    // Draw Y-axis
    glColor3ub(0, 255, 0);
    // 線分を描画します。2024/03 CJS
    renderer.clear();
    renderer.appendLine(-0.1 * scl, -length * scl * 0.5, -0.1 * scl, -0.1 * scl, length * scl, -0.1 * scl);
    renderer.drawLines();

    glPushMatrix();
    arrows[1] = gluNewQuadric();
    gluQuadricDrawStyle(arrows[1], GLU_FILL);
    glTranslated(-0.1 * scl, length * scl, -0.1 * scl);
    glRotated(-90.0f, 1, 0, 0);
    gluCylinder(arrows[1], length / 30 * scl, 0.0f, length / 5 * scl, 8, 8);
    gluDeleteQuadric(arrows[1]);
    glPopMatrix();

    // Draw Z-axis
    glColor3ub(0, 0, 255);
    // 線分を描画します。2024/03 CJS
    renderer.clear();
    renderer.appendLine(-0.1 * scl, -0.1 * scl, -length * scl * 0.5, -0.1 * scl, -0.1 * scl, length * scl);
    renderer.drawLines();

    glPushMatrix();
    arrows[2] = gluNewQuadric();
    gluQuadricDrawStyle(arrows[2], GLU_FILL);
    glTranslated(-0.1 * scl, -0.1 * scl, length * scl);
    gluCylinder(arrows[2], length / 30 * scl, 0.0f, length / 5 * scl, 8, 8);
    gluDeleteQuadric(arrows[2]);
    glPopMatrix();

    glColor3ub(255, 255, 255);
    glDisable(GL_DEPTH_TEST);
    //glPrint3d(length*1.1*scl, 0.0, 0.0, (void *)font, "X");
    //glPrint3d(0.0, length*1.1*scl, 0.0, (void *)font, "Y");
    //glPrint3d(0.0, 0.0, length*1.1*scl, (void *)font, "Z");
    glEnable(GL_DEPTH_TEST);
}

void cnt_switch_visibility_by_wall()
{
    double dsmin = 0.2e-10, dsmax = 2.0e-10;
    double dsminsq = dsmin * dsmin;
    double dsmaxsq = dsmax * dsmax;
    bool debug = true;
    double* crx, * cry; // coordinates of points for cross sectional shape
    crx = new double[atom.natom + 1];
    cry = new double[atom.natom + 1];
    bool* flg; // flag for search
    flg = new bool[atom.natom + 1];
    double rr, dss; // for calculation
    for (int i = 0; i <= atom.natom; i++) { flg[i] = true; }
    bool* flg2;
    flg2 = new bool[atom.natom + 1];
    for (int i = 0; i <= atom.natom; i++) { flg2[i] = true; }

    // calculate center (x,y)
    double center_x, center_y;
    calc_center_xy(center_x, center_y);

    for (int iwall = 0; iwall < 6; iwall++) { // iwall=0 is outermost wall

        // remove crowded points (flg = false means out of search list)
        for (int i = 1; i <= atom.natom - 1; i++) {
            if (flg[i] && flg2[i]) {
                for (int j = i + 1; j <= atom.natom; j++) {
                    if (flg[j] && flg2[j]) { // if i & j are too close, j is removed from search list
                        dss = dsquare(atom.rx[j] - atom.rx[i]) + dsquare(atom.ry[j] - atom.ry[i]);
                        if (dss < dsminsq) { flg[j] = false; }
                    }
                }
            }
        }
        // count number of points
        int npts = 0;
        for (int i = 1; i <= atom.natom; i++) { if (flg[i] && flg2[i]) npts++; }
        // find farthest atom (which should be in outermost wall)
        int ifar = 0; double rrfar = 0;
        for (int i = 1; i <= atom.natom; i++) {
            if (flg[i] && flg2[i]) {
                rr = dsquare(atom.rx[i] - center_x) + dsquare(atom.ry[i] - center_y);
                if (rr > rrfar) { rrfar = rr; ifar = i; }
            }
        }

        // (for debugging) output line shape
        /*
        if (debug) {
          FILE *fp = fopen("cnt_cross_section_shape.d","w");
          for (int i=1; i<atom.natom; i++) {
            if (flg[i]) fprintf(fp,"%f %f\n",atom.rx[i]*1.0e10,atom.ry[i]*1.0e10);
          }
          fclose(fp);
        }
        */

        // find points for outermost wall (starting from ifar)
        flg[ifar] = false;
        int itag = ifar; int npts_out = 1; int inn; double xnn, ynn;
        crx[1] = atom.rx[ifar]; cry[1] = atom.ry[ifar];
        for (int k = 1; k <= npts; k++) { // k is dummy
            double rrmin = 1.0e10;
            for (int i = 1; i <= atom.natom; i++) {
                if (flg[i] && flg2[i]) {
                    dss = dsquare(atom.rx[i] - atom.rx[itag]) + dsquare(atom.ry[i] - atom.ry[itag]);
                    if (dss < rrmin) { rrmin = dss; inn = i; xnn = atom.rx[i]; ynn = atom.ry[i]; }
                }
            }
            if (rrmin < dsmaxsq) {
                flg[inn] = false; npts_out++; crx[npts_out] = xnn; cry[npts_out] = ynn; itag = inn;
            }
            else {
                break;
            }
        }
        crx[0] = crx[npts_out]; cry[0] = cry[npts_out];

        // switch visibility for atoms in outermost wall
        for (int i = 1; i <= atom.natom; i++) {
            if (flg2[i]) { atom.visible[i] = false; }
            for (int j = 1; j <= npts_out; j++) {
                rr = dsquare(atom.rx[i] - crx[j]) + dsquare(atom.ry[i] - cry[j]);
                if (rr < dsmaxsq) { // if yes, atom is in iwall-th outermost wall
                    //atom.visible[i] = true;
                    if (show_cnt_layer[iwall]) { atom.visible[i] = true; }
                    flg2[i] = false; // exclude atom from list
                }
            }
        }

        // (for debugging) output line shape
        /*
        if (debug) {
          FILE *fpp = fopen("cnt_cross_section_shape_outermost_sort.d","w");
          for (int i=0; i<=npts_out; i++) {
            fprintf(fpp,"%f %f\n",crx[i]*1.0e10,cry[i]*1.0e10);
          }
          fclose(fpp);
        }
        */

    }

    delete[] crx; delete[] cry; delete[] flg; delete[] flg2;
}
