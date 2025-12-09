#pragma once

#include "BaseWindow.h"

/**
 * Set Parameter ウィンドウです。
 **/
class SetParameterWindow : public BaseWindow
{
public:
    /* コンストラクターです。*/
    SetParameterWindow();

    /**
     * イベント ハンドラーです。
     * @param event イベント識別子
     **/
    int handle(int event) override;
};
