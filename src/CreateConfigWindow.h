#pragma once

#include <memory>
#include "BaseWindow.h"
#include "EditAtomWindow.h"

/**
 * 設定ファイルを作成するウィンドウです。
 **/
class CreateConfigWindow : public BaseWindow
{
public:
    /* コンストラクターです。*/
    CreateConfigWindow();

    /**
     * イベント ハンドラーです。
     * @param event イベント識別子
     **/
    int handle(int event) override;

private:
    /* Edit atom ウィンドウです。*/
    std::unique_ptr<EditAtomWindow> edit_atom_window;
};
