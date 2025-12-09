#pragma once

#include <FL/Fl.H>
#include <FL/Fl_Choice.H>
#include <string>
#include <vector>

/**
 * ドロップダウン リスト ウィジェットです。
 * ウィジェットの値を登録した外部変数と紐づけ、
 * 外部変数の値の変化がウィジェットに反映されます。
 **/
class CustomChoice : public Fl_Choice {
public:

    /**
     * コンストラクターです。
     * @param X ウィジェットの x 座標
     * @param Y ウィジェットの y 座標
     * @param W ウィジェットの幅
     * @param H ウィジェットの高さ
     * @param L ウィジェットのラベル
     **/
    CustomChoice(int X, int Y, int W, int H, const char* L = 0);
    /* デストラクターです。*/
    virtual ~CustomChoice();

    /* 入力と変数を紐づけます。*/
    void attach(int* selected_index);
    /* 変数の更新を開始します。*/
    void start_polling();
    /* 変数の更新を終了します。*/
    void stop_polling();
    /* 変数の更新を行うコールバック関数です。*/
    static void polling_callback(void* v);

private:
    /* 選択されたインデックスを表す外部変数のポインターです。*/
    int* _external_selected_index = nullptr;
    /* 変数の更新を行っているかどうかです。*/
    bool _polling = false;
    /* 変数の同期を行います。*/
    void sync_with_external_variables();

    /**
     * イベント ハンドラーです。
     * @param event イベント識別子
     **/
    int handle(int event);
};
