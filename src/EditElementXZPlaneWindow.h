#pragma once

#include "BaseWindow.h"

/**
 * Edit element XZ plane ウィンドウです。
 **/
class EditElementXZPlaneWindow : public BaseWindow
{
public:
    /* コンストラクターです。*/
    EditElementXZPlaneWindow();

    /**
     * イベント ハンドラー
     * @param event イベント識別子
     **/
    int handle(int event) override;
};
