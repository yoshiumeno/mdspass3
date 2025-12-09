#define _USE_MATH_DEFINES
#include <cmath>



#include "MeshRenderer.h"

/**
 * ベクトル (0, 0, 1) をベクトル (x, y, z) に回転させる回転行列を求めます。
 * @param x                 [in] ベクトルの X 成分です。
 * @param y                 [in] ベクトルの Y 成分です。
 * @param z                 [in] ベクトルの Z 成分です。
 * @param rotation_matrix   [out] 回転行列を受け取るサイズ 9 の配列です。
 * */
inline void getRotationMatrix(float x, float y, float z, float rotation_matrix[9])
{
    // 回転軸ベクトルを求めます。
    float kx = -y;
    float ky = x;
    float kz = 0;

    // 回転角を計算します。
    float length = sqrt(x * x + y * y + z * z);
    float cos_theta = z / length;
    float sin_theta = sqrt(1 - cos_theta * cos_theta);

    // 回転行列を計算します。
    rotation_matrix[0] = cos_theta + kx * kx * (1 - cos_theta);
    rotation_matrix[1] = kx * ky * (1 - cos_theta) - kz * sin_theta;
    rotation_matrix[2] = kx * kz * (1 - cos_theta) + ky * sin_theta;
    rotation_matrix[3] = ky * kx * (1 - cos_theta) + kz * sin_theta;
    rotation_matrix[4] = cos_theta + ky * ky * (1 - cos_theta);
    rotation_matrix[5] = ky * kz * (1 - cos_theta) - kx * sin_theta;
    rotation_matrix[6] = kz * kx * (1 - cos_theta) - ky * sin_theta;
    rotation_matrix[7] = kz * ky * (1 - cos_theta) + kx * sin_theta;
    rotation_matrix[8] = cos_theta + kz * kz * (1 - cos_theta);
}

/**
 * ベクトルに回転行列を適用します。
 * @param rotation_matrix   [in] 回転行列です。
 * @param x0                [in] 回転対象のベクトルの X 成分です。
 * @param y0                [in] 回転対象のベクトルの Y 成分です。
 * @param z0                [in] 回転対象のベクトルの Z 成分です。
 * @param x1                [out] 回転後のベクトルの X 成分です。
 * @param y1                [out] 回転後のベクトルの Y 成分です。
 * @param z1                [out] 回転後のベクトルの Z 成分です。
 * */
inline void applyRotationMatrix(float rotation_matrix[9], float x0, float y0, float z0, float& x1, float& y1, float& z1)
{
    x1 = rotation_matrix[0] * x0 + rotation_matrix[1] * y0 + rotation_matrix[2] * z0;
    y1 = rotation_matrix[3] * x0 + rotation_matrix[4] * y0 + rotation_matrix[5] * z0;
    z1 = rotation_matrix[6] * x0 + rotation_matrix[7] * y0 + rotation_matrix[8] * z0;
}


// メッシュ、線分の登録を全て解除します。
void MeshRenderer::clear()
{
    mesh_vertex_count = 0;
    mesh_face_count = 0;
    line_vertex_count = 0;
    line_count = 0;
}

// 球体のメッシュを描画対象として追加します。
void MeshRenderer::appendSphere(
    GLfloat x,
    GLfloat y,
    GLfloat z,
    GLfloat r,
    const GLfloat* color,
    int division)
{
    // 緯度と経度方向のメッシュの分割数を決定します。
    const int num_vertices_per_sphere = (division + 1) * division;
    const int num_faces_per_sphere = division * division * 2;

    // 球体メッシュの頂点数を計算し、記録先の頂点配列のサイズが足りない場合は拡張します。
    const int new_vertice_count = mesh_vertex_count + num_vertices_per_sphere;
    if (mesh_vertex_list.size() < new_vertice_count) mesh_vertex_list.resize(new_vertice_count);
    if (mesh_color_list.size() < new_vertice_count)  mesh_color_list.resize(new_vertice_count);
    if (mesh_normal_list.size() < new_vertice_count) mesh_normal_list.resize(new_vertice_count);

    // 頂点の座標、色、法線ベクトルを配列に登録します。
    for (GLuint i = 0; i <= division; i++) {
        const float theta = i * M_PI / division;
        for (GLuint j = 0; j < division; j++) {
            // 角度を決定します。
            const float phi = j * 2 * M_PI / division;
            // 配列上の登録先インデックスを決定します。
            const int index = mesh_vertex_count + i * division + j;

            // 頂点の座標を登録します。
            mesh_vertex_list[index].x = r * sin(theta) * cos(phi) + x;
            mesh_vertex_list[index].y = r * sin(theta) * sin(phi) + y;
            mesh_vertex_list[index].z = r * cos(theta) + z;
            // 頂点の色を登録します。
            mesh_color_list[index].r = color[0];
            mesh_color_list[index].g = color[1];
            mesh_color_list[index].b = color[2];
            mesh_color_list[index].a = color[3];
            // 頂点の法線ベクトルを登録します。
            mesh_normal_list[index].x = sin(theta) * cos(phi);
            mesh_normal_list[index].y = sin(theta) * sin(phi);
            mesh_normal_list[index].z = cos(theta);
        }
    }

    // 球体メッシュの面数を計算し、記録先の面データ配列のサイズが足りない場合は拡張します。
    const int new_face_count = mesh_face_count + num_faces_per_sphere;
    if (mesh_face_list.size() < new_face_count) mesh_face_list.resize(new_face_count);

    // 面データを配列に登録します。
    for (GLuint i = 0; i < division; i++) {
        for (GLuint j = 0; j < division; j++) {
            // 配列上の登録先インデックスを決定します。
            const int base_index = mesh_face_count + i * division * 2 + j * 2;

            // 四角形を対角線で分割した最初の三角形を構成する頂点のインデックスを登録します。
            mesh_face_list[base_index].v1 = mesh_vertex_count + i * division + j;
            mesh_face_list[base_index].v2 = mesh_vertex_count + (i + 1) * division + j;
            mesh_face_list[base_index].v3 = mesh_vertex_count + (i + 1) * division + (j + 1 == division ? 0 : j + 1);

            // 四角形を対角線で分割した二番目の三角形を構成する頂点のインデックスを登録します。
            mesh_face_list[base_index + 1].v1 = mesh_vertex_count + i * division + j;
            mesh_face_list[base_index + 1].v2 = mesh_vertex_count + (i + 1) * division + (j + 1 == division ? 0 : j + 1);
            mesh_face_list[base_index + 1].v3 = mesh_vertex_count + i * division + (j + 1 == division ? 0 : j + 1);
        }
    }

    // 頂点数と面数を管理する変数を更新します。
    mesh_vertex_count = new_vertice_count;
    mesh_face_count = new_face_count;
}

// 矢印のメッシュを描画対象として追加します。
void MeshRenderer::appendArrowVector(
    GLfloat x, GLfloat y, GLfloat z,
    GLfloat vx, GLfloat vy, GLfloat vz,
    float scale)
{
    // 矢印始点 (円錐の反対側) を scale で補正します。
    x *= scale;
    y *= scale;
    z *= scale;

    // ベクトル長を scale で補正します。
    vx *= scale;
    vy *= scale;
    vz *= scale;

    // 矢印の方向ベクトルと長さを計算します。
    GLfloat length = sqrt(vx * vx + vy * vy + vz * vz);
    GLfloat ux = vx / length;
    GLfloat uy = vy / length;
    GLfloat uz = vz / length;

    // 矢印の底面の回転に使用する回転行列を求めます。
    float rotmat[9];
    getRotationMatrix(ux, uy, uz, rotmat);

    // 矢印の軸部分と円錐部分の長さを決定します。
    const GLfloat line_length = length * ARROW_SHAFT_LENGTH_RATIO;
    const GLfloat cone_length = length * ARROW_CONE_LENGTH_RATIO;

    // 矢印メッシュの頂点数を計算し、記録先の頂点配列のサイズが足りない場合は拡張します。
    const int vertex_count = ARROW_SLICES * 3 + 2;
    if (mesh_vertex_list.size() < mesh_vertex_count + vertex_count) {
        mesh_vertex_list.resize(mesh_vertex_count + vertex_count);
    }

    // 矢印の軸部分の頂点座標を配列に登録します。
    for (int i = 0; i < ARROW_SLICES; i++)
    {
        const int shaft_base_index = mesh_vertex_count + i;

        // 矢印の中心軸に対する、頂点の相対座標を計算します。
        GLfloat angle = 2 * M_PI * i / ARROW_SLICES;
        GLfloat rx = ARROW_SHAFT_RADIUS * cos(angle);
        GLfloat ry = ARROW_SHAFT_RADIUS * sin(angle);

        float px, py, pz;
        applyRotationMatrix(rotmat, rx, ry, 0, px, py, pz);

        // 軸の頂点座標を配列に登録します。
        mesh_vertex_list[shaft_base_index].x = x + px;
        mesh_vertex_list[shaft_base_index].y = y + py;
        mesh_vertex_list[shaft_base_index].z = z + pz;

        mesh_vertex_list[shaft_base_index + ARROW_SLICES].x = x + px + ux * line_length;
        mesh_vertex_list[shaft_base_index + ARROW_SLICES].y = y + py + uy * line_length;
        mesh_vertex_list[shaft_base_index + ARROW_SLICES].z = z + pz + uz * line_length;
    }

    // 矢印の円錐部分の頂点座標を配列に登録します。
    for (int i = 0; i < ARROW_SLICES; i++)
    {
        const int cone_base_index = mesh_vertex_count + i + ARROW_SLICES * 2;

        // 円錐の軸に対する、頂点の相対座標を計算します。
        GLfloat angle = 2 * M_PI * i / ARROW_SLICES;
        GLfloat rx = ARROW_CONE_RADIUS * cos(angle);
        GLfloat ry = ARROW_CONE_RADIUS * sin(angle);

        float px, py, pz;
        applyRotationMatrix(rotmat, rx, ry, 0, px, py, pz);

        // 円錐の頂点座標を配列に登録します。
        mesh_vertex_list[cone_base_index].x = x + px + ux * line_length;
        mesh_vertex_list[cone_base_index].y = y + py + uy * line_length;
        mesh_vertex_list[cone_base_index].z = z + pz + uz * line_length;
    }

    // 円錐の先端は全て同じ座標なので、頂点配列には 1 つだけ登録します。
    mesh_vertex_list[mesh_vertex_count + ARROW_SLICES * 3].x = x + vx;
    mesh_vertex_list[mesh_vertex_count + ARROW_SLICES * 3].y = y + vy;
    mesh_vertex_list[mesh_vertex_count + ARROW_SLICES * 3].z = z + vz;


    // 矢印メッシュの面数を計算し、記録先の面データの配列のサイズが足りない場合は拡張します。
    const int face_count = ARROW_SLICES * 2 + ARROW_SLICES * 1;
    if (mesh_face_list.size() < mesh_face_count + face_count) {
        mesh_face_list.resize(mesh_face_count + face_count);
    }

    // 矢印の軸部分の面を登録します。
    for (int i = 0; i < ARROW_SLICES; i++)
    {
        const int base_index = mesh_face_count + i;
        const int i1 = (i + 1 == ARROW_SLICES) ? 0 : i + 1;

        mesh_face_list[base_index].v1 = mesh_vertex_count + i;
        mesh_face_list[base_index].v2 = mesh_vertex_count + i1;
        mesh_face_list[base_index].v3 = mesh_vertex_count + ARROW_SLICES + i;

        mesh_face_list[base_index + ARROW_SLICES].v1 = mesh_vertex_count + i1;
        mesh_face_list[base_index + ARROW_SLICES].v2 = mesh_vertex_count + ARROW_SLICES + i1;
        mesh_face_list[base_index + ARROW_SLICES].v3 = mesh_vertex_count + ARROW_SLICES + i;
    }

    // 矢印の円錐部分の面を登録します。
    for (int i = 0; i < ARROW_SLICES; i++)
    {
        const int cone_Base_index = mesh_face_count + i + ARROW_SLICES * 2;
        const int i1 = (i + 1 == ARROW_SLICES) ? 0 : i + 1;

        mesh_face_list[cone_Base_index].v1 = mesh_vertex_count + ARROW_SLICES * 3;
        mesh_face_list[cone_Base_index].v2 = mesh_vertex_count + ARROW_SLICES * 2 + i;
        mesh_face_list[cone_Base_index].v3 = mesh_vertex_count + ARROW_SLICES * 2 + i1;
    }

    // 頂点数と面数を管理する変数を更新します。
    mesh_vertex_count += vertex_count;
    mesh_face_count += face_count;
}

// 矢印のメッシュを描画対象として追加します。
void MeshRenderer::appendArrowPosition(
    GLfloat x0, GLfloat y0, GLfloat z0,
    GLfloat x1, GLfloat y1, GLfloat z1)
{
    appendArrowVector(x0, y0, z0, x1 - x0, y1 - y0, z1 - z0);
}

// 三角形を描画対象として追加します。
void MeshRenderer::appendTriangle(
    GLfloat x0, GLfloat y0, GLfloat z0,
    GLfloat x1, GLfloat y1, GLfloat z1,
    GLfloat x2, GLfloat y2, GLfloat z2)
{
    // 頂点配列と面データ配列の必要なサイズを計算し、足りなければ各配列を拡張します。
    const int vertex_count = 3;
    if (mesh_vertex_list.size() < mesh_vertex_count + vertex_count) {
        mesh_vertex_list.resize(mesh_vertex_count + vertex_count);
    }
    const int face_count = 1;
    if (mesh_face_list.size() < mesh_face_count + face_count) {
        mesh_face_list.resize(mesh_face_count + face_count);
    }

    // 頂点を登録します。
    mesh_vertex_list[mesh_vertex_count].x = x0;
    mesh_vertex_list[mesh_vertex_count].y = y0;
    mesh_vertex_list[mesh_vertex_count].z = z0;
    mesh_vertex_list[mesh_vertex_count + 1].x = x1;
    mesh_vertex_list[mesh_vertex_count + 1].y = y1;
    mesh_vertex_list[mesh_vertex_count + 1].z = z1;
    mesh_vertex_list[mesh_vertex_count + 2].x = x2;
    mesh_vertex_list[mesh_vertex_count + 2].y = y2;
    mesh_vertex_list[mesh_vertex_count + 2].z = z2;

    // 面を登録します。
    mesh_face_list[mesh_face_count].v1 = mesh_vertex_count;
    mesh_face_list[mesh_face_count].v2 = mesh_vertex_count + 1;
    mesh_face_list[mesh_face_count].v3 = mesh_vertex_count + 2;

    // 配列上の要素数を出力用変数に登録します。
    mesh_vertex_count += vertex_count;
    mesh_face_count += face_count;
}

// 線分を 1 つ配列に追加します。
void MeshRenderer::appendLine(
    GLfloat x0, GLfloat y0, GLfloat z0,
    GLfloat x1, GLfloat y1, GLfloat z1)
{
    // 頂点配列と線分配列の必要なサイズを計算し、各配列のサイズが足りない場合は拡張します。
    const int vertex_count = 2;
    if (line_vertex_list.size() < line_vertex_count + vertex_count) {
        line_vertex_list.resize(line_vertex_count + vertex_count);
    }
    const int _line_count = 1;
    if (line_list.size() < line_count + _line_count) {
        line_list.resize(line_count + _line_count);
    }

    // 頂点を登録します。
    line_vertex_list[line_vertex_count].x = x0;
    line_vertex_list[line_vertex_count].y = y0;
    line_vertex_list[line_vertex_count].z = z0;
    line_vertex_list[line_vertex_count + 1].x = x1;
    line_vertex_list[line_vertex_count + 1].y = y1;
    line_vertex_list[line_vertex_count + 1].z = z1;

    // 線を登録します。
    line_list[line_count].v1 = line_vertex_count;
    line_list[line_count].v2 = line_vertex_count + 1;

    // 頂点数と面数を管理する変数を更新します。
    line_vertex_count += vertex_count;
    line_count += _line_count;
}

// 円柱のメッシュを配列に追加します。
void MeshRenderer::appendCylinder(
    GLfloat x0,
    GLfloat y0,
    GLfloat z0,
    GLfloat x1,
    GLfloat y1,
    GLfloat z1,
    float radius,
    int div_number)
{
    // 太さをオリジナルの処理に合わせて 1/100 (定数 CYLINDER_RADIUS_CORRECT_RATIO で定義) にします。
    radius *= CYLINDER_RADIUS_CORRECT_RATIO;

    // 方向ベクトルを計算します。
    // 長さは 0.9 倍 (定数 ARROW_SHAFT_RATIO で定義) しておきます。
    GLfloat vx = (x1 - x0) * ARROW_SHAFT_LENGTH_RATIO;
    GLfloat vy = (y1 - y0) * ARROW_SHAFT_LENGTH_RATIO;
    GLfloat vz = (z1 - z0) * ARROW_SHAFT_LENGTH_RATIO;
    GLfloat length = sqrt(vx * vx + vy * vy + vz * vz);
    GLfloat ux = vx / length;
    GLfloat uy = vy / length;
    GLfloat uz = vz / length;

    // 円柱の回転に使用する回転行列を求めます。
    float rotmat[9];
    getRotationMatrix(ux, uy, uz, rotmat);

    // 円柱の頂点数を計算し、記録先の頂点配列のサイズが足りない場合は拡張します。
    const int vertex_count = (div_number) * 2;
    if (mesh_vertex_list.size() < mesh_vertex_count + vertex_count) {
        mesh_vertex_list.resize(mesh_vertex_count + vertex_count);
    }
    if (mesh_normal_list.size() < mesh_vertex_count + vertex_count) {
        mesh_normal_list.resize(mesh_vertex_count + vertex_count);
    }

    // 頂点座標と法線ベクトルを配列に登録します。
    for (int i = 0; i < div_number; i++)
    {
        const int base_index = mesh_vertex_count + i;

        // 円柱の軸に対する、頂点の相対座標を計算します。
        GLfloat angle = 2 * M_PI * i / div_number;
        GLfloat rx = cos(angle);
        GLfloat ry = sin(angle);

        float px, py, pz;
        applyRotationMatrix(rotmat, rx, ry, 0, px, py, pz);

        // 円柱の頂点座標と法線ベクトルを配列に登録します。
        mesh_vertex_list[base_index].x = x0 + px * radius;
        mesh_vertex_list[base_index].y = y0 + py * radius;
        mesh_vertex_list[base_index].z = z0 + pz * radius;
        mesh_normal_list[base_index].x = px;
        mesh_normal_list[base_index].y = py;
        mesh_normal_list[base_index].z = pz;

        mesh_vertex_list[base_index + div_number].x = x0 + vx + px * radius;
        mesh_vertex_list[base_index + div_number].y = y0 + vy + py * radius;
        mesh_vertex_list[base_index + div_number].z = z0 + vz + pz * radius;
        mesh_normal_list[base_index + div_number].x = px;
        mesh_normal_list[base_index + div_number].y = py;
        mesh_normal_list[base_index + div_number].z = pz;
    }


    // 面数を計算し、記録先の面データ配列のサイズが足りない場合は拡張します。
    const int face_count = div_number * 2;
    if (mesh_face_list.size() < mesh_face_count + face_count) {
        mesh_face_list.resize(mesh_face_count + face_count);
    }

    // 面を登録します。
    for (int i = 0; i < div_number; i++)
    {
        const int base_face_index = mesh_face_count + i;
        int i1 = (i + 1 == div_number ? 0 : i + 1);

        mesh_face_list[base_face_index].v1 = mesh_vertex_count + i;
        mesh_face_list[base_face_index].v2 = mesh_vertex_count + i1;
        mesh_face_list[base_face_index].v3 = mesh_vertex_count + div_number + i;

        mesh_face_list[base_face_index + div_number].v1 = mesh_vertex_count + i1;
        mesh_face_list[base_face_index + div_number].v2 = mesh_vertex_count + div_number + i1;
        mesh_face_list[base_face_index + div_number].v3 = mesh_vertex_count + div_number + i;
    }

    // 頂点数と面数を管理する変数を更新します。
    mesh_vertex_count += vertex_count;
    mesh_face_count += face_count;
}


// 配列に登録されているメッシュを描画します。
void MeshRenderer::drawMeshesWithColor()
{
    // ライティングと深度テスト、カリングを有効にします。
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    // 頂点毎の色情報を有効にします。
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);

    // 使用可能な配列の種類を指定します。
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);

    // 描画に使用する配列を指定します。
    glVertexPointer(3, GL_FLOAT, 0, mesh_vertex_list.data());
    glColorPointer(4, GL_FLOAT, 0, mesh_color_list.data());
    glNormalPointer(GL_FLOAT, sizeof(Normal), mesh_normal_list.data());

    // メッシュを描画します。
    glDrawElements(GL_TRIANGLES, mesh_face_count * 3, GL_UNSIGNED_INT, mesh_face_list.data());

    // この描画固有の設定を無効にします。
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
}

// 配列に登録されているメッシュを描画します。
void MeshRenderer::drawMeshesWithNormal()
{
    // 深度テスト、カリングを有効にします。
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    // 使用可能な配列の種類を指定します。
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);

    // 描画に使用する配列を指定します。
    glVertexPointer(3, GL_FLOAT, 0, mesh_vertex_list.data());
    glNormalPointer(GL_FLOAT, sizeof(Normal), mesh_normal_list.data());

    // メッシュを描画します。
    glDrawElements(GL_TRIANGLES, mesh_face_count * 3, GL_UNSIGNED_INT, mesh_face_list.data());

    // この描画固有の設定を無効にします。
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
}

// 配列に登録されているメッシュを描画します。
void MeshRenderer::drawMeshes()
{
    // 深度テスト、カリングを有効にします。
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    // 使用可能な配列の種類を指定します。
    glEnableClientState(GL_VERTEX_ARRAY);

    // 描画に使用する配列を指定します。
    glVertexPointer(3, GL_FLOAT, 0, mesh_vertex_list.data());

    // メッシュを描画します。
    glDrawElements(GL_TRIANGLES, mesh_face_count * 3, GL_UNSIGNED_INT, mesh_face_list.data());

    // この描画固有の設定を無効にします。
    glDisableClientState(GL_VERTEX_ARRAY);
}

// 配列に登録されているメッシュをワイヤーフレームとして描画します。
void MeshRenderer::drawWireframe()
{
    // ライティングと深度テストを有効にします。
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);

    // ポリゴンの塗りつぶしモードを線に設定します。
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    // 頂点毎の色情報を有効にします。
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);

    // 使用可能な配列の種類を指定します。
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);

    // 頂点配列と法線配列にデータを設定します。
    glVertexPointer(3, GL_FLOAT, sizeof(Vertex), mesh_vertex_list.data());
    glNormalPointer(GL_FLOAT, sizeof(Normal), mesh_normal_list.data());

    // 面を描画するために、glDrawElements 関数を呼び出します。
    glDrawElements(GL_TRIANGLES, mesh_face_count * 3, GL_UNSIGNED_INT, mesh_face_list.data());

    // この描画固有の設定を無効にします。
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisable(GL_COLOR_MATERIAL);

    // ポリゴンの塗りつぶしモードを元に戻します。
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

// 配列に登録されている線分を描画します。
void MeshRenderer::drawLines(GLfloat line_width)
{
    // 線の太さを指定します。
    glLineWidth(line_width);

    // 深度テストのみ有効にします。
    glEnable(GL_DEPTH_TEST);

    // 使用可能な配列の種類を指定します。
    glEnableClientState(GL_VERTEX_ARRAY);

    // 描画に使用する配列を指定します。
    glVertexPointer(3, GL_FLOAT, 0, line_vertex_list.data());

    // 線分を描画します。
    glDrawElements(GL_LINES, line_count * 2, GL_UNSIGNED_INT, line_list.data());

    // 頂点配列を無効にします。
    glDisableClientState(GL_VERTEX_ARRAY);

    // 線の太さを標準に戻します。
    glLineWidth(1.0);
}
