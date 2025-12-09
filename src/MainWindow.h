#pragma once

#include <memory>
#include "BaseWindow.h"
#include "SetParameterWindow.h"
#include "StressWindow.h"
#include "ExtraWindow.h"
#include "EditElementWindow.h"
#include "EditElementXZPlaneWindow.h"
#include "CreateConfigWindow.h"
#include "PotentialArgWindow.h"
#include "MdViewerWindow.h"
#include "OpenFileWindow.h"
#include "QuitWindow.h"

/**
 * メイン ウィンドウです。
 **/
class MainWindow : public BaseWindow
{
public:
    /* コンストラクターです。*/
    MainWindow();

    /**
     * イベント ハンドラーです。
     * @param event イベント識別子
     **/
    int handle(int event) override;

private:
    /* Set Parameter ウィンドウです。*/
    std::unique_ptr<SetParameterWindow> set_parameter_window;
    /* Stress ウィンドウです。*/
    std::unique_ptr<StressWindow> stress_window;
    /* Edit element ウィンドウです。*/
    std::unique_ptr<EditElementWindow> edit_element_window;
    /* Edit element xz plane ウィンドウです。*/
    std::unique_ptr<EditElementXZPlaneWindow> edit_element_xz_plane_window;
    /* Extra ウィンドウです。*/
    std::unique_ptr<ExtraWindow> extra_window;
    /* Create config ウィンドウです。*/
    std::unique_ptr<CreateConfigWindow> create_config_window;
    /* Potential arg ウィンドウです。*/
    std::unique_ptr<PotentialArgWindow> potential_arg_window;
    /* MD viewer ウィンドウです。*/
    std::unique_ptr<MdViewerWindow> md_viewer_window;
    /* Open File ウィンドウです。*/
    std::unique_ptr<OpenFileWindow> open_file_window;
    // Quit ウィンドウ追加 (YU2025.01)
    std::unique_ptr<QuitWindow> quit_window;
};
