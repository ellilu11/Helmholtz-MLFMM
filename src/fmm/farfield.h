#pragma once

#include "node.h"

class FMM::Farfield {

public:
    Farfield(const std::shared_ptr<Node>&);

    void buildMpoleCoeffs(const std::shared_ptr<Node>&);

    void mergeMpoleCoeffs(const std::shared_ptr<Node>&);

    void buildLocalCoeffs(const std::shared_ptr<Node>&);

    void evaluateSols();

private:
    void buildLevels();

    void buildRadPats();

    void buildRadPatsOfNode(const std::shared_ptr<Node>&);

    void resizeCoeffs(const std::shared_ptr<Node>&);

    Coeffs getShiftedLocalCoeffs(Node*, int) const;

    void translateCoeffs(const std::shared_ptr<Node>&);

    void evalFarSols(const std::shared_ptr<Node>&);

    static void addInterpCoeffs(const Coeffs&, Coeffs&, int, int);

    static void addAnterpCoeffs(const Coeffs&, Coeffs&, int, int);
};