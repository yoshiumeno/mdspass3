#include <FL/Fl.H>
#if defined(_WIN32) || defined(_WIN64)
#include <Windows.h>
#endif
#include <FL/gl.h>
#include <FL/glu.h>
#include <FL/glut.H>
#include "GlWindow.h"
#include "myglutdisplay.h"
#include "myheader.h"

using namespace std;

// OpenGL の描画に関する関数のプロトタイプ宣言です。
void writedata_initialize();
void view_save();
void view_read();
void rot2qua(float q[], float r[]);
void qmul(float r[], const float p[], const float q[]);
void qua2rot(float r[], float q[]);
void set_atom_color();
void set_atom_color(int ia);

//extern void ensure_display_lists(); // Dec 2025

// bugfix for Mac (Dec 2025)
// OpenGLのコンテキストがcurrentになっていない時にglDeleteLists/glGenListsなどを
// 呼んでいる箇所があり（readconfig.cppなど）、Win/Linuxでは「たまたま」問題に
// ならなかったが、Macではsegmentation faultが出てしまう。この不具合を修正。
static int  g_num_display_lists = 0; // 現在確保しているリスト数
static bool g_display_lists_dirty = true; // 「作り直しが必要」フラグ
// readconfig から呼ぶ用：「次回描画時に display list を作り直して」
void mark_display_lists_dirty()
{
    g_display_lists_dirty = true;
}
// GlWindow から呼ぶ用：必要ならここで display list を再生成
void ensure_display_lists()
{
    if (!g_display_lists_dirty)
       return;
    // コンテキストはGlWindows側で必ずcurrentになっている前提
    if (objects != 0 && g_num_display_lists > 0 && glIsList(objects)) {
        glDeleteLists(objects, g_num_display_lists);
    }
    g_num_display_lists = atom.natom * 3;
    if (g_num_display_lists > 0) {
        objects = glGenLists(g_num_display_lists);
    } else {
        objects = 0;
    }
    g_display_lists_dirty = false;
}

GlWindow::GlWindow(int x_, int y_, int w_, int h_, const char* l)
    : Fl_Gl_Window(x_, y_, w_, h_, l)
{
    // OpenGL の初期化です。
    // バッファ モードを設定します。
    mode(FL_RGB | FL_DOUBLE | FL_DEPTH);
//    mode(FL_RGB | FL_SINGLE | FL_DEPTH);

    // Dec 2025: glEnable moved to draw()
    //glEnable(GL_LIGHTING);
    //glEnable(GL_LIGHT0);
    //glEnable(GL_DEPTH_TEST);
}

// マウスやキーボードからの入力を処理する関数です。
int GlWindow::handle(int event)
{
    static GLuint selection[SELECTIONS];
    static GLint hits = 0;

    switch (event)
    {
        // マウス操作に関する処理です。
        // マウス ボタンが押されたイベントです。
        case FL_PUSH:
            // ウィンドウにフォーカスし、キーボード入力を処理できるようにします。
            Fl::focus(this);

            // 現在のマウスの座標を取得します。
            cx = Fl::event_x();
            cy = Fl::event_y();

            // 右ボタンが押された場合です。
            if (Fl::event_button() == FL_RIGHT_MOUSE)
            {
                // 現在のスケール値を保存します。
                cs = scl;
                // オブジェクトの現在位置を保存します。
                cpos0 = obj_pos[0];
                cpos1 = obj_pos[1];
                // 回転行列からクォータニオンに変換します。
                rot2qua(cq, rotmat);
                // 右ボタンが押されたフラグを立てます。
                right_button = true;
            }
            // 左ボタンが押された場合です。
            else if (Fl::event_button() == FL_LEFT_MOUSE)
            {
                // ピッキング処理を行います。
                perform_picking(cx, cy, selection, hits);
            }
            return 1;
        // マウス ボタンが離されたイベントです。
        case FL_RELEASE:
            // 右ボタンの場合です。
            if (Fl::event_button() == FL_RIGHT_MOUSE)
            {
                // クォータニオンを更新します。
                cq[0] = tq[0];
                cq[1] = tq[1];
                cq[2] = tq[2];
                cq[3] = tq[3];
                // 右ボタンが離されたフラグを立てます。
                right_button = false;
            }
            // 左ボタンの場合です。
            else if (Fl::event_button() == FL_LEFT_MOUSE)
            {
                // 左ボタンのリリース処理を行います。
                end_picking(cx, cy, selection, hits);
            }
            return 1;
        // マウス ドラッグ イベントです。
        case FL_DRAG:
        {
            int x = Fl::event_x();
            int y = Fl::event_y();
            if (right_button)
            {
                if (Fl::event_state(FL_CTRL))
                { // Ctrl キーが押されている場合です。
                    double dx = ( x - cx ) * sx;
                    double dy = ( y - cy ) * sy;
                    scl = cs * ( 1 - dy );
                }
                else if (Fl::event_state(FL_SHIFT))
                { // Shift キーが押されている場合です。
                    double dx = ( x - cx ) * sx;
                    double dy = ( y - cy ) * sy;
                    obj_pos[0] = cpos0 + dx;
                    obj_pos[1] = cpos1 - dy;
                }
                else
                { // それ以外の場合 (Ctrl キーや Shift キーが押されていない) です。
                    float dx = ( x - cx ) * sx;
                    float dy = ( y - cy ) * sy;
                    float a = sqrtf(dx * dx + dy * dy);
                    if (a != 0.0)
                    {
                        float ar = a * 2.0 * M_PI * 0.5;
                        float as = sin(ar) / a;
                        float dq[4] = { cos(ar), dy * as, dx * as, 0.0 };
                        qmul(tq, dq, cq);
                        qua2rot(rotmat, tq); // クォータニンから回転行列に変換します。
                    }
                }
                // 再描画します。
                // Timer イベントで自動更新されるため、コメントアウトしました。
                //redraw();
            }
            return 1;
        }
        // キーボード操作に関する処理です。
        case FL_KEYDOWN:
        {
            // 入力されたキーに応じて処理を行います。
            int key = Fl::event_key();

            if (key == 'e')
            {
                // Cycle through ensemble types
                ++ensemble %= MAXENSTYPE;
                return 1;
            }
            else if (key == 'r')
            {
                org_x = 0.0; org_y = 0.0; org_z = 0.0;
                rotationX = 0.0; rotationY = 0.0;
                for (int i = 0; i <= 15; i++) { rotmat[i] = 0; }
                rotmat[0] = 1; rotmat[6] = -1; rotmat[9] = 1; rotmat[15] = 1;
                obj_pos[0] = 0; obj_pos[1] = 0; obj_pos[2] = 0; scl = 1.0;
                writedata_initialize();

            }
            else if (key == 'V')
            {
                view_save();
            }
            else if (key == 'v')
            {
                view_read();
            }
            return 1;
        }
        default:
            return Fl_Gl_Window::handle(event);
    };
    // Timer イベントで自動更新されるため、コメントアウトしました。
    //redraw();
}

void GlWindow::draw()
{
    if (!valid())
    {
        resize(x(), y(), w(), h());

        // Dec 2025
        change_ortho(x(), y(), w(), h());
        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);
        glEnable(GL_DEPTH_TEST);

    }

    // Dec 2025
    // display list が汚れていればここで作り直す（コンテキストは必ず current）
    ensure_display_lists();

    myGlutDisplay();

    // ダブル バッファー モードで描画します。
    // NOTE: 2024/08 CJS
    // 実機インストールした Ubuntu 22.04 で正常に動作しなかったため、コメントアウトしています。
    //swap_buffers();
}

void GlWindow::change_ortho(int x, int y, int w, int h)
{
    // 2024/02 CJS main.cpp myGlutReshape 関数から移動
    xy_aspect = (float)w / (float)h;
    glViewport(0, 0, w, h);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if (ortho_state)
    {
        GLfloat light0_ambient[] = { 0.1f, 0.1f, 0.3f, 1.0f };
        GLfloat light0_diffuse[] = { .6f, .6f, 1.0f, 1.0f };
        GLfloat light0_position[] = { 200.0f, 200.0f, 200 - 6.0f, 1.0f };
        glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
        glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
        glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
        glOrtho(-w / 200.0, w / 200.0, -h / 200.0, h / 200.0, -100.0, 100.0);
    }
    else
    {
        GLfloat light0_ambient[] = { 0.1f, 0.1f, 0.3f, 1.0f };
        GLfloat light0_diffuse[] = { .6f, .6f, 1.0f, 1.0f };
        GLfloat light0_position[] = { 20.0f, 20.0f, 200 - 4.0f, 1.0f };
        glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
        glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
        glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
        gluPerspective(30.0, xy_aspect, 1.0, 100.0);
        glTranslated(0.0, 0.0, -1.0);
    }
    glMatrixMode(GL_MODELVIEW);

    // 2024/02 CJS main.cpp myGlutReshape 関数から移動 ここまで
}

void GlWindow::resize(int x, int y, int w, int h)
{
    Fl_Gl_Window::resize(x, y, w, h);

    // 遠近法あり、なしの描画モードを切り替えます。
    // Dec 2025: moved to draw()
    //change_ortho(x, y, w, h);

    size_w = w; size_h = h;

    // マウス ドラッグの感度を設定します。
    // 2024/02 CJS マウス感度を上げました。
    // 元の値は以下です。 
    // sx = 1.0 / (double)w 
    // sy = 1.0 / (double)h
    sx = 1.5 / (double)w;
    sy = 1.5 / (double)h;

    // キャプチャ用の幅と、高さを更新します。
    capture_width = w;
    capture_height = h;
}

void GlWindow::perform_picking(int x, int y, GLuint* selection, GLint hits)
{
    // 2024/02 CJS main.cpp myGlutMouse 関数から移動
    GLuint* ptr; GLint vp[4];
    glSelectBuffer(SELECTIONS, selection);
    glDisable(GL_LIGHTING);//YU20200520 (bugfix for Ubuntu20.04)
    glRenderMode(GL_SELECT);
    glMatrixMode(GL_PROJECTION);
    glInitNames();
    glPushName(-1);
    glPushMatrix();
    glLoadIdentity();
    glGetIntegerv(GL_VIEWPORT, vp);
    gluPickMatrix(x, vp[3] - y - 1, 1, 1, vp);
    ensure_display_lists(); // Dec 2025
    if (ortho_state)
    {
        glOrtho(-size_w / 200.0, size_w / 200.0, -size_h / 200.0, size_h / 200.0, -100.0, 100.0);
    }
    else
    {
        gluPerspective(30.0, (double)vp[2] / (double)vp[3], 1.0, 100.0);
        glTranslated(0.0, 0.0, -1.0);
    }
    glMatrixMode(GL_MODELVIEW);
    for (int i = 1; i <= atom.natom; i++)
    {
        glLoadName(i);
        glCallList(objects + i);
    }
    if (draw_replica > 0)
    {
        for (int i = 1; i <= icnt; i++)
        {
            glLoadName(ibase + i);
            glCallList(ibase + i);
        }
    }
    glEnable(GL_LIGHTING);//YU20200520
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    hits = glRenderMode(GL_RENDER);
    ptr = selection;

    //printf("Atom information: hits = %d\n", hits);
    if (hits > 0)
    {
        unsigned int j, n = ptr[0];
        ptr += 3;
        for (j = 0; j < n; j++)
        {
            if (*ptr <= atom.natom)
            {
                atom.Calcq(*ptr);
                printf("Atom = %d [%s] Pos (A): %f %f %f\n",
                    *ptr, atom.asp[*ptr], atom.rx[*ptr] / ang, atom.ry[*ptr] / ang, atom.rz[*ptr] / ang);
                printf("                (frac): %f %f %f\n",
                    atom.qx[*ptr], atom.qy[*ptr], atom.qz[*ptr]);
                printf("               F (eV/A): %f %f %f (total)\n",
                    atom.fx[*ptr] / eV * ang, atom.fy[*ptr] / eV * ang, atom.fz[*ptr] / eV * ang);
                printf("                       : %f %f %f (from loading)\n",
                    atom.fx_l[*ptr] / eV * ang, atom.fy_l[*ptr] / eV * ang, atom.fz_l[*ptr] / eV * ang);
                printf("               Ene (eV): %f ,  Mass (kg): %e\n",
                    atom.epot[*ptr] / eV, atom.wm[*ptr]);
                printf("               V  (m/s): %f %f %f\n",
                    atom.vx[*ptr], atom.vy[*ptr], atom.vz[*ptr]);
                printf("           Stress (MPa): xx = %f yy = %f zz = %f\n",
                    atom.satom[*ptr][0][0] / 1e6, atom.satom[*ptr][1][1] / 1e6,
                    atom.satom[*ptr][2][2] / 1e6);
                printf("                         xy = %f yz = %f zx = %f\n",
                    atom.satom[*ptr][0][1] / 1e6, atom.satom[*ptr][1][2] / 1e6,
                    atom.satom[*ptr][2][0] / 1e6);
                printf("    Number of neighbors: %d\n", book.alistnum[*ptr]);
                if (createconfig_mode||editatom_mode)
                {
                    iatom_pick = *ptr;
                    strcpy(config_atom, atom.asp[iatom_pick]);
                    atomrx = atom.rx[iatom_pick] / ang;
                    atomry = atom.ry[iatom_pick] / ang;
                    atomrz = atom.rz[iatom_pick] / ang;
                }
            }
            else
            {
                int ia = iatom[*ptr - atom.natom - objects];
                printf("Replica Atom = %d  Pos: %f %f %f  Rep-index: %d\n",
                    ia, atom.rx[ia] * 1e9, atom.ry[ia] * 1e9, atom.rz[ia] * 1e9,
                    repidx[*ptr - objects]);
            }

            if (edit_elem_mode > 0)
            {
                int iimax = 4; if (edit_elem_mode == 2) { iimax = 3; }
                int ifnd = 0;
                for (int ii = 0; ii < iimax; ii++)
                {
                    int ix = *ptr; int iy = 0;
                    if (ix > atom.natom) { ix = iatom[*ptr - atom.natom - objects]; iy = repidx[*ptr - atom.natom - objects]; }
                    if (( select_atom[ii] == ix ) && ( select_atom_repidx[ii] == iy ))
                    {
                        if (*ptr <= atom.natom)
                        {
                            memcpy(::color[*ptr], yellow, sizeof(GLfloat) * 4);
                        }
                        else
                        {
                            memcpy(::color[*ptr - objects], yellow, sizeof(GLfloat) * 4);
                        }
                        select_atom[ii] = 0; select_atom_repidx[ii] = 0;
                        ifnd = 1;
                        break;
                    }
                }
                if (ifnd == 0)
                {
                    for (int ii = 0; ii < iimax; ii++)
                    {
                        if (select_atom[ii] == 0)
                        {
                            if (*ptr <= atom.natom)
                            {
                                memcpy(::color[*ptr], green, sizeof(GLfloat) * 4);
                            }
                            else
                            {
                                memcpy(::color[*ptr - objects], green, sizeof(GLfloat) * 4);
                            }
                            if (*ptr <= atom.natom)
                            {
                                select_atom[ii] = *ptr; select_atom_repidx[ii] = 0;
                            }
                            else
                            {
                                select_atom[ii] = iatom[*ptr - atom.natom - objects];
                                select_atom_repidx[ii] = repidx[*ptr - atom.natom - objects];
                            }
                            printf("Selected Atom %d = %d \n", ii, select_atom[ii]);
                            break;
                        }
                    }
                }
            }
            else if (measure_mode > 0)
            {
                int iimax = 3;
                int ifnd = 0;
                for (int ii = 0; ii < iimax; ii++)
                {
                    int ix = *ptr; int iy = 0;
                    if (ix > atom.natom)
                    {
                        ix = iatom[*ptr - atom.natom - objects];
                        iy = repidx[*ptr - atom.natom - objects];
                    }
                    if (( select_atom[ii] == ix ) && ( select_atom_repidx[ii] == iy ))
                    {
                        if (*ptr <= atom.natom)
                        {
                            set_atom_color(*ptr);
                        }
                        else { set_atom_color(*ptr - objects); }

                        select_atom[ii] = 0; select_atom_repidx[ii] = 0;
                        ifnd = 1; break;
                    }
                }
                if (ifnd == 0)
                {
                    for (int ii = 0; ii < iimax; ii++)
                    {
                        if (select_atom[ii] == 0)
                        {
                            if (*ptr <= atom.natom)
                            {
                                memcpy(::color[*ptr], green, sizeof(GLfloat) * 4);
                            }
                            else
                            {
                                memcpy(::color[*ptr - objects], green, sizeof(GLfloat) * 4);
                            }
                            if (*ptr <= atom.natom)
                            {
                                select_atom[ii] = *ptr; select_atom_repidx[ii] = 0;
                            }
                            else
                            {
                                select_atom[ii] = iatom[*ptr - atom.natom - objects];
                                select_atom_repidx[ii] = repidx[*ptr - atom.natom - objects];
                            }
                            printf("Selected Atom %d = %d \n", ii, select_atom[ii]);
                            break;
                        }
                    }
                }

            }
            else
            { // if not edit_elem_mode nor measure_mode
                memcpy(::color[*ptr], blue, sizeof(GLfloat) * 4);
            } // end if edit_elem_mode & measure_mode
        }
    }
    // 2024/02 CJS main.cpp myGlutMouse 関数から移動 ここまで

    // Timer イベントで自動更新されるため、コメントアウトしました。 
    //GlWindow::redraw();
}

void GlWindow::end_picking(int x, int y, GLuint* selection, GLint hits)
{
    // 2024/02 CJS main.cpp myGlutMouse 関数から移動
    GLuint* ptr = selection;
    GLint vp[4];
    glSelectBuffer(SELECTIONS, selection);
    glDisable(GL_LIGHTING);//YU20200520 (bugfix for Ubuntu20.04)
    glRenderMode(GL_SELECT);
    glMatrixMode(GL_PROJECTION);
    glInitNames();
    glPushName(-1);
    glPushMatrix();
    glLoadIdentity();
    glGetIntegerv(GL_VIEWPORT, vp);
    gluPickMatrix(x, vp[3] - y - 1, 1, 1, vp);
    ensure_display_lists(); // Dec 2025
    if (ortho_state)
    {
        glOrtho(-size_w / 200.0, size_w / 200.0, -size_h / 200.0, size_h / 200.0, -100.0, 100.0);
    }
    else
    {
        gluPerspective(30.0, (double)vp[2] / (double)vp[3], 1.0, 100.0);
        glTranslated(0.0, 0.0, -1.0);
    }
    glMatrixMode(GL_MODELVIEW);
    for (int i = 1; i <= atom.natom; i++)
    {
        glLoadName(i);
        glCallList(objects + i);
    }
    if (draw_replica > 0)
    {
        for (int i = 1; i <= icnt; i++)
        {
            glLoadName(ibase + i);
            glCallList(ibase + i);
        }
    }
    glEnable(GL_LIGHTING);//YU20200520
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    hits = glRenderMode(GL_RENDER);

    if (hits > 0)
    {
        unsigned int j, n = ptr[0];
        ptr += 3;
        for (j = 0; j < n; j++)
        {
            if (( edit_elem_mode > 0 ) || ( measure_mode > 0 ))
            {

            }
            else
            {   // return to original color
                memcpy(::color[*ptr], color0[*ptr], sizeof(GLfloat) * 4);
            }
        }
    }
    // 2024/02 CJS main.cpp myGlutMouse 関数から移動 ここまで
}
